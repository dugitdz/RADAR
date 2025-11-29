clc; clear; close all;
%% ===================== PARÂMETROS CONFIGURÁVEIS ===================== %%
N        = 16;
overlap  = 0.8;
hop      = max(1, round(N*(1-overlap)));
win_smooth_raw   = 128;
win_smooth_final = 128;
f_low  = 0.8;
f_high = 3.0;

% === PARÂMETROS DE TRACKING ===
bpm_change_limit = 10;  % Mudança máxima esperada entre janelas
penalty_factor   = 2.5; % Penalidade para saltos bruscos

%% ========================== LEITURA ========================== %%
data  = readmatrix('C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phase.csv');
tpuro  = data(:,1);
heart  = data(:,4);

%% ====== TEMPO ======
t  = (tpuro - tpuro(1)) / 1000;
dt = mean(diff(t));
srate = 1/dt;
nt = numel(t);

%% ====== SEGMENTAÇÃO POR BURACOS ======
dt_all   = diff(t);
break_ix = find(dt_all > 0.5);
if isempty(break_ix)
    seg_starts = 1;
    seg_ends   = nt;
else
    seg_starts = [1; break_ix + 1];
    seg_ends   = [break_ix; nt];
end
nSeg = numel(seg_starts);

%% ====== FREX PELO ESPECTRO DO HEART (ORIGINAL) ======
Nfft = 2*nt - 1;
f_fft = linspace(0, srate, Nfft);
H = fft(heart, Nfft);
[~, imax] = max(abs(H));
frex = f_fft(imax);

%% ====== DEFINIR h = X / frex (ORIGINAL) ======
X = 0.6;               
h = X / frex;
h = max(h, 0.1);
h = min(h, 5.0);

%% ========================== WAVELET & CONVOLUÇÃO (ORIGINAL) ========================== %%
idx   = (1:nt) - ceil(nt/2);
t_aux = idx(:) * dt;
gaus = exp((-4*log(2)*t_aux.^2)/h^2) .* exp(1i*2*pi*frex*t_aux);
nGaus  = length(gaus);
nConv  = nGaus + nt - 1;
heartX = fft(heart, nConv);
gausX  = fft(gaus,  nConv);
gausX  = gausX ./ max(gausX);
convX     = heartX .* gausX;
convs_all = ifft(convX);
startIdx = floor(nGaus/2) + 1;
endIdx   = startIdx + nt - 1;
convs    = convs_all(startIdx:endIdx);

%% ====== FFT POR JANELAS COM TRACKING ======
halfN   = N/2;
rel_idx = -halfN : (halfN-1);

t_window      = [];
freq_h_window = [];
k = 0;

for s = 1:nSeg
    seg_start = seg_starts(s);
    seg_end   = seg_ends(s);
    center_list = seg_start : hop : seg_end;
    
    last_valid_bpm = NaN; % Reset do tracking por segmento

    for center_theoretical = center_list
        idx_theo = center_theoretical + rel_idx;
        idx_win  = idx_theo(idx_theo >= seg_start & idx_theo <= seg_end);
        
        if numel(idx_win) < 4, continue; end
        
        t_seg = t(idx_win);
        x_seg = real(convs(idx_win));
        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end
        fs_seg = 1/mean(dt_seg);
        x_seg = x_seg - mean(x_seg);
        
        % --- APLICAÇÃO DE TAPER (HANNING) CORRIGIDA ---
        win_taper = hanning(numel(x_seg));
        x_seg_tapered = x_seg(:) .* win_taper; % (:) garante vetor coluna
        
        Xw   = fft(x_seg_tapered);
        Lwin = numel(x_seg_tapered);
        half = floor(Lwin/2);
        f_ax = (0:half-1) * fs_seg / Lwin;
        
        magX = abs(Xw(1:half));
        mask = (f_ax >= f_low) & (f_ax <= f_high);
        
        if ~any(mask), continue; end
        
        % Candidatos
        f_cands = f_ax(mask);
        m_cands = magX(mask);
        bpm_cands = f_cands * 60;
        
        % --- LÓGICA DE TRACKING ---
        if isnan(last_valid_bpm)
            % Sem memória: Pega o pico mais forte
            [~, imax_win] = max(m_cands);
            if length(imax_win) > 1, imax_win = imax_win(1); end
            bpm_chosen = bpm_cands(imax_win);
        else
            % Com memória: Score Ponderado
            dist_bpm = abs(bpm_cands - last_valid_bpm);
            weights  = 1 ./ (1 + penalty_factor * (dist_bpm / bpm_change_limit).^2);
            score    = m_cands .* weights;
            
            [~, best_idx] = max(score);
            if length(best_idx) > 1, best_idx = best_idx(1); end
            bpm_chosen = bpm_cands(best_idx);
        end
        
        last_valid_bpm = bpm_chosen;
        
        k = k + 1;
        t_window(k,1)      = mean(t_seg);
        freq_h_window(k,1) = bpm_chosen;
    end
end

%% ====== 1) INTERPOLAÇÃO ======
freq_h_raw = nan(nt,1);
idx_good_win = ~isnan(freq_h_window);

if any(idx_good_win)
    freq_h_raw = interp1(t_window(idx_good_win), freq_h_window(idx_good_win), ...
                         t, 'pchip', nan); 
end

%% ====== 2) SUAVIZAÇÃO POR SEGMENTO ======
freq_h_smooth = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_h_smooth(idx_seg) = movmean(freq_h_raw(idx_seg), win_smooth_raw, 'omitnan');
end

%% ====== 3) PREENCHIMENTO DE GAPS ======
idx_good = ~isnan(freq_h_smooth);
if any(idx_good)
    t_targets = t(~idx_good);
    freq_h_smooth(find(~idx_good)) = interp1(t(idx_good), freq_h_smooth(idx_good), ...
                                             t_targets, 'linear', 'extrap');
end

%% ====== 4) SUAVIZAÇÃO FINAL ======
freq_h_final = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_h_final(idx_seg) = movmean(freq_h_smooth(idx_seg), win_smooth_final, 'omitnan');
end

%% ====== GRÁFICO FINAL ======
figure;
plot(t_window, freq_h_window, 'o', 'Color', [0.7 0.7 0.7]); hold on;
plot(t, freq_h_raw, 'b.', 'MarkerSize', 4);
plot(t, freq_h_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('HR (bpm)');
title(sprintf('HR Final — Wavelet + Tracking (h=%.3f)', h));
grid on;
legend('Tracking','Raw','Final');

%% ====== SALVAR RESULTADO ======
freq_h       = freq_h_raw;
freq_h_cont  = freq_h_final;
mean_h_cont  = freq_h_final;
freq_br      = nan(size(t));
freq_br_cont = freq_br;
mean_br_cont = freq_br;

save('freq_results.mat', 't', ...
    'freq_h_raw','freq_h', ...
    'freq_h_cont','mean_h_cont', ...
    'freq_br','freq_br_cont','mean_br_cont', ...
    't_window','freq_h_window', ...
    'seg_starts','seg_ends', ...
    'N','overlap','hop', ...
    'win_smooth_raw','win_smooth_final', ...
    'frex','h','X');