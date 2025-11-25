clc; clear; close all;

%% ========================== PATH ==========================
phases_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv';

%% ========================== LEITURA ==========================
data   = readmatrix(phases_path,'CommentStyle', '#', 'Delimiter', ',');
tpuro  = data(:,1);      % tempo bruto (ms)
total  = data(:,2);
breath = data(:,3);
heart  = data(:,4);

%% ========================== NORMALIZAÇÃO DE FASE ==========================
breath = mod(breath + pi, 2*pi) - pi;
heart  = mod(heart  + pi, 2*pi) - pi;

% tempo em segundos começando em zero
t  = (tpuro - tpuro(1)) / 1000;
nt = numel(t);

%% ========================== CHECAR GAPS EM t ==========================
dt = diff(t);
gap_idx = find(dt > 0.5);   % gaps maiores que 0.5 s

if ~isempty(gap_idx)
    fprintf('Foram encontrados %d gaps com dt > 0.5 s:\n', numel(gap_idx));
    for k = 1:numel(gap_idx)
        i = gap_idx(k);
        fprintf('  Entre i = %d (t = %.3f s) e i = %d (t = %.3f s): dt = %.3f s\n', ...
            i, t(i), i+1, t(i+1), dt(i));
    end
else
    fprintf('Nenhum dt > 0.5 s encontrado.\n');
end

%% ========================== PARÂMETROS FFT (JANELA PADRÃO) ==========================
N       = 64;        % tamanho da janela
Nmin    = N/2;       % tamanho mínimo aceitável
overlap = 0.9;       % sobreposição (90%)
hop     = max(1, round(N*(1 - overlap)));   % passo em amostras

if nt < Nmin
    error('Poucos pontos para fazer FFT com janela (nt < Nmin).');
end

%% ========================== ACUMULADORES ==========================
sum_freq_br = zeros(size(t));
sum_w_br    = zeros(size(t));
sum_freq_h  = zeros(size(t));
sum_w_h     = zeros(size(t));

t_window       = [];
freq_h_window  = [];
freq_br_window = [];

%% ========================== LOOP FFT DESLIZANTE ==========================
g_start   = 1;
janela_id = 0;

while (g_start + N - 1) <= nt
    g_end = g_start + N - 1;
    idx   = g_start:g_end;
    Lraw  = numel(idx);

    if Lraw < Nmin
        break;
    end

    % janela atual
    t_win  = t(idx);
    br_win = breath(idx);
    h_win  = heart(idx);

    dt_win = diff(t_win);
    if any(dt_win <= 0) || any(isnan(dt_win))
        % janela ruim, só pula
        g_start = g_start + hop;
        continue;
    end

    fs = 1 / mean(dt_win);

    % ==================== FFT ====================
    Br       = fft(br_win);
    H        = fft(h_win);
    nFFT     = numel(br_win);
    half_fft = floor(nFFT/2);
    f_axis   = (0:half_fft-1) * fs / nFFT;

    [~,kbr] = max(abs(Br(1:half_fft)));
    [~,kh ] = max(abs(H(1:half_fft)));

    fbr_hz = f_axis(kbr);
    fh_hz  = f_axis(kh);

    % salvar ponto central da janela em bpm
    t_window(end+1,1)       = mean(t_win);
    freq_br_window(end+1,1) = fbr_hz * 60;
    freq_h_window(end+1,1)  = fh_hz  * 60;

    % acumular frequência em todos os pontos da janela
    sum_freq_br(idx) = sum_freq_br(idx) + fbr_hz;
    sum_w_br(idx)    = sum_w_br(idx)    + 1;
    sum_freq_h(idx)  = sum_freq_h(idx)  + fh_hz;
    sum_w_h(idx)     = sum_w_h(idx)     + 1;

    % ===== PRINT SIMPLES DA JANELA =====
    janela_id = janela_id + 1;
    fprintf('Janela %3d: t_ini = %.3f s | t_fim = %.3f s | N = %d\n', ...
        janela_id, t_win(1), t_win(end), numel(idx));

    % próximo início
    g_start = g_start + hop;
end

%% ========================== FREQ BRUTA POR AMOSTRA ==========================
freq_br = nan(size(t));
freq_h  = nan(size(t));

mask_br = sum_w_br > 0;
mask_h  = sum_w_h  > 0;

freq_br(mask_br) = sum_freq_br(mask_br) ./ sum_w_br(mask_br);
freq_h(mask_h)   = sum_freq_h(mask_h)   ./ sum_w_h(mask_h);

% Hz -> bpm
freq_br = freq_br * 60;
freq_h  = freq_h  * 60;

%% ========================== CURVA CONTÍNUA GLOBAL ==========================
valid_h  = ~isnan(freq_h);
valid_br = ~isnan(freq_br);

freq_h_cont  = interp1(t(valid_h),  freq_h(valid_h),  t, 'pchip', 'extrap');
freq_br_cont = interp1(t(valid_br), freq_br(valid_br), t, 'pchip', 'extrap');

%% ========================== MOVMEAN ==========================
win_smooth = 256;

mean_h_cont  = movmean(freq_h_cont,  win_smooth);
mean_br_cont = movmean(freq_br_cont, win_smooth);

%% ========================== PLOTS ==========================
figure;
plot(t, freq_h_cont, '-', 'LineWidth', 1.0); hold on;
plot(t, mean_h_cont, 'r', 'LineWidth', 1.5);
xlabel('Tempo (s)'); ylabel('HR (bpm)');
title('HR contínua + movmean'); grid on;

figure;
plot(t, freq_br_cont, '-', 'LineWidth', 1.0); hold on;
plot(t, mean_br_cont, 'r', 'LineWidth', 1.5);
xlabel('Tempo (s)'); ylabel('RR (bpm)');
title('RR contínua + movmean'); grid on;

%% ========================== SAVE ==========================
save('freq_results.mat', 't', ...
    'freq_h', 'freq_br', ...
    'freq_h_cont', 'freq_br_cont', ...
    'mean_h_cont', 'mean_br_cont');

save('freq_windows.mat','t_window','freq_h_window','freq_br_window');
