clc; clear; close all;

%% ========================== PATH ==========================
phases_path = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phase.csv';

%% ========================== LEITURA ==========================
data   = readmatrix(phases_path,'CommentStyle', '#', 'Delimiter', ',');
tpuro  = data(:,1);      % tempo bruto (ms)
total  = data(:,2);
breath = data(:,3);
heart  = data(:,4);
heart = unwrap(heart);
%% ========================== NORMALIZAÇÃO DE FASE ==========================
% limite a fase para [-pi, +pi] sem alterar FFT, só para corrigir wraps muito grandes
breath = mod(breath + pi, 2*pi) - pi;
heart  = mod(heart  + pi, 2*pi) - pi;

% tempo em segundos começando em zero
t  = (tpuro - tpuro(1)) / 1000;
nt = numel(t);
srate = 1/mean(diff(t));

%% ========================== PARÂMETROS FFT ==========================
N       = 128;
Nmin    = N/4;
overlap = 0.99;
hop     = max(1, round(N*(1-overlap)));

gap_thr = 0.5;

% filtros
br_low  = 0.10;
br_high = 0.50;
h_low   = 0.80;
h_high  = 2.50;

%% ========================== SEGMENTAÇÃO ==========================
dt      = diff(t);
gap_idx = find(dt > gap_thr);

seg_start = [1; gap_idx+1];
seg_end   = [gap_idx; nt];
numSeg    = numel(seg_start);

%% ========================== ACUMULADORES ==========================
sum_freq_br = zeros(size(t));
sum_w_br    = zeros(size(t));
sum_freq_h  = zeros(size(t));
sum_w_h     = zeros(size(t));

t_window       = [];
freq_h_window  = [];
freq_br_window = [];

%% ========================== LOOP DE FFT ==========================
for s = 1:numSeg
    
    i1 = seg_start(s);
    i2 = seg_end(s);
    L  = i2 - i1 + 1;
    if L < Nmin, continue; end
    
    first_end = i1 + Nmin - 1;
    if first_end > i2, continue; end
    
    last_start = [];
    
    for g_end = first_end : hop : i2
        
        g_start = max(i1, g_end - (N-1));
        idx_raw = g_start:g_end;
        Lraw    = numel(idx_raw);
        if Lraw < Nmin, continue; end
        
        p = floor(log2(Lraw));
        M = 2^p;
        if M < Nmin, continue; end
        
        if Lraw == M
            start2 = g_start;
            end2   = g_end;
        else
            mid = round((g_start + g_end)/2);
            start2 = mid - floor((M-1)/2);
            end2   = start2 + M - 1;
            if start2 < i1, start2 = i1; end
            if end2 > i2, end2 = i2; start2 = end2 - M + 1; end
        end
        
        if start2 < i1 || end2 > i2 || (end2-start2+1)~=M
            continue;
        end
        
        if ~isempty(last_start) && start2 <= last_start
            hopM = max(1, round(M*(1-overlap)));
            new_start2 = last_start + hopM;
            new_end2   = new_start2 + M - 1;
            if new_end2 > i2, continue; end
            start2 = new_start2;
            end2   = new_end2;
        end
        
        last_start = start2;
        
        idx    = start2:end2;
        t_win  = t(idx);
        br_raw = breath(idx);
        h_raw  = heart(idx);
        
        dt_win = diff(t_win);
        if any(dt_win <= 0) || any(isnan(dt_win)), continue; end
        
        fs = 1/mean(dt_win);
        
        % filtros
        if fs > 2*br_high
            [b,a] = butter(4,[br_low br_high]/(fs/2));
            br_win = filtfilt(b,a,br_raw);
        else
            br_win = br_raw;
        end
        
        if fs > 2*h_high
            [b,a] = butter(4,[h_low h_high]/(fs/2));
            h_win = filtfilt(b,a,h_raw);
        else
            h_win = h_raw;
        end
        
        % -------- FFT --------
        % %%% VERSÃO ANTIGA (CORTAVA METADE)
        % Br       = fft(br_win);
        % H        = fft(h_win);
        % nFFT     = numel(br_win);
        % half_fft = floor(nFFT/2);
        % f_axis   = (0:half_fft-1) * fs / nFFT;
        % [~,kbr] = max(abs(Br(1:half_fft)));
        % [~,kh ] = max(abs(H(1:half_fft)));
        % fbr_hz = f_axis(kbr);
        % fh_hz  = f_axis(kh);
        
        % %%% ALTERADO AQUI: usa espectro completo, sem cortar "frequência negativa"
        Br   = fft(br_win);
        H    = fft(h_win);
        nFFT = numel(br_win);
        f_axis = (0:nFFT-1) * fs / nFFT;   % eixo 0..(fs-df)

        [~,kbr] = max(abs(Br));            % pico em todo o espectro
        [~,kh ] = max(abs(H));

        fbr_hz = f_axis(kbr);
        fh_hz  = f_axis(kh);
        
        % salvar
        t_window(end+1,1)       = mean(t_win);
        freq_br_window(end+1,1) = fbr_hz*60;
        freq_h_window(end+1,1)  = fh_hz*60;
        
        % acumular
        sum_freq_br(idx) = sum_freq_br(idx) + fbr_hz;
        sum_w_br(idx)    = sum_w_br(idx) + 1;
        sum_freq_h(idx)  = sum_freq_h(idx) + fh_hz;
        sum_w_h(idx)     = sum_w_h(idx) + 1;
    end
end

%% ========================== FREQ BRUTA POR AMOSTRA ==========================
freq_br = nan(size(t));
freq_h  = nan(size(t));

mask_br = sum_w_br > 0;
mask_h  = sum_w_h  > 0;

freq_br(mask_br) = sum_freq_br(mask_br)./sum_w_br(mask_br);
freq_h(mask_h)   = sum_freq_h(mask_h)./sum_w_h(mask_h);

freq_br = freq_br*60;
freq_h  = freq_h*60;

%% ========================== CURVA CONTÍNUA GLOBAL ==========================
valid_h  = ~isnan(freq_h);
valid_br = ~isnan(freq_br);

freq_h_cont  = interp1(t(valid_h),  freq_h(valid_h),  t, 'pchip', 'extrap');
freq_br_cont = interp1(t(valid_br), freq_br(valid_br), t, 'pchip', 'extrap');

%% ========================== SUAVIZAÇÃO GAUSSIANA (O QUE VOCÊ ESTAVA FAZENDO) ==========================
h = 2;
n = length(t);
t_deslocado = t - t(floor((n+1)/2));

gaus = exp(-4*log(2)*t_deslocado.^2/h^2);
gaus = gaus/sum(gaus);

ngaus   = length(gaus);
nsignal = length(freq_h_cont);
nConv   = ngaus + nsignal - 1;
half    = floor(ngaus/2);

G = fft(gaus,       nConv);
S = fft(freq_h_cont,nConv);
freq_h_new = ifft(G .* S);

freq_h_new = freq_h_new(half+1:end-half);

% ================= AJUSTE PEDIDO =================
% Se for menor que 60, arredonda para 60 bpm
freq_h_new(freq_h_new < 60) = 60;
  % mantém como você colocou

%% ========================== MOVMEAN ==========================
win_smooth = 256;

mean_h_cont  = movmean(freq_h_cont,  win_smooth);
mean_br_cont = movmean(freq_br_cont, win_smooth);

%% ========================== PLOTS ESSENCIAIS ==========================
figure;
plot(t, freq_h_cont, '-', 'LineWidth', 1.0); hold on;
plot(t, freq_h_new, 'm');
plot(t, mean_h_cont, 'r', 'LineWidth', 1.5);
xlabel('Tempo (s)'); ylabel('HR (bpm)');
title('HR contínua + gauss + movmean');
grid on;

figure;
plot(t, freq_br_cont, '-', 'LineWidth', 1.0); hold on;
plot(t, mean_br_cont, 'r', 'LineWidth', 1.5);
xlabel('Tempo (s)'); ylabel('RR (bpm)');
title('RR contínua + movmean');
grid on;
%mean_h_cont = freq_h_new;
%% ========================== SAVE ==========================
save('freq_results.mat', 't', ...
    'freq_h', 'freq_br', ...
    'freq_h_cont', 'freq_br_cont', ...
    'mean_h_cont', 'mean_br_cont');

save('freq_windows.mat','t_window','freq_h_window','freq_br_window');
