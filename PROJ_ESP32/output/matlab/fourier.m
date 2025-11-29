clc; clear; close all;

%% ===================== PARÂMETROS CONFIGURÁVEIS ===================== %%
N        = 16;
overlap  = 0.8;
hop      = max(1, round(N*(1-overlap)));

win_smooth_raw   = 128;
win_smooth_final = 128;

% banda fisiológica do coração
f_low_h  = 0.8;
f_high_h = 3.0;

% banda fisiológica da respiração (ajusta se quiser)
f_low_br  = 0.1;
f_high_br = 0.7;

%% ========================== LEITURA ========================== %%
data   = readmatrix('C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phases.csv');
tpuro  = data(:,1);
breath = data(:,3);
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

%% ====== FREX PELO ESPECTRO DO HEART ======
Nfft  = 2*nt - 1;
f_fft = linspace(0, srate, Nfft);
H     = fft(heart, Nfft);
[~, imax_h] = max(abs(H));
frex_h = f_fft(imax_h);

%% ====== DEFINIR h_h = X / frex_h ======
X = 0.6;               % <<<<—— FIXO AQUI
h_h = X / frex_h;
h_h = max(h_h, 0.1);
h_h = min(h_h, 5.0);

%% ========================== WAVELET & CONVOLUÇÃO (HEART) ========================== %%
idx   = (1:nt) - ceil(nt/2);
t_aux = idx(:) * dt;

gaus_h = exp((-4*log(2)*t_aux.^2)/h_h^2) .* exp(1i*2*pi*frex_h*t_aux);

nGaus_h  = length(gaus_h);
nConv_h  = nGaus_h + nt - 1;

heartX = fft(heart, nConv_h);
gausX_h  = fft(gaus_h,  nConv_h);
gausX_h  = gausX_h ./ max(gausX_h);

convX_h     = heartX .* gausX_h;
convs_all_h = ifft(convX_h);

startIdx_h = floor(nGaus_h/2) + 1;
endIdx_h   = startIdx_h + nt - 1;
convs_h    = convs_all_h(startIdx_h:endIdx_h);

%% ====== FFT POR JANELAS (HEART, SEGMENTADO) ======
halfN   = N/2;
rel_idx = -halfN : (halfN-1);

freq_weighted_sum_h = zeros(nt,1);
weight_sum_h        = zeros(nt,1);

t_window_h      = [];
freq_h_window   = [];
k_h = 0;

for s = 1:nSeg
    seg_start = seg_starts(s);
    seg_end   = seg_ends(s);

    center_list = seg_start : hop : seg_end;

    for center_theoretical = center_list
        idx_theo = center_theoretical + rel_idx;
        idx_win  = idx_theo(idx_theo >= seg_start & idx_theo <= seg_end);

        if numel(idx_win) < 4
            continue;
        end

        t_seg = t(idx_win);
        x_seg = real(convs_h(idx_win));

        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end

        fs_seg = 1/mean(dt_seg);

        x_seg = x_seg - mean(x_seg);

        Xw   = fft(x_seg);
        Lwin = numel(x_seg);
        half = floor(Lwin/2);

        f_ax = (0:half-1) * fs_seg / Lwin;
        magX = abs(Xw(1:half));

        mask = (f_ax >= f_low_h) & (f_ax <= f_high_h);
        if ~any(mask), continue; end

        [~, imax_win] = max(magX(mask));
        f_est_vec = f_ax(mask);
        f_est_h   = f_est_vec(imax_win);

        bpm_h = f_est_h * 60;

        k_h = k_h + 1;
        t_window_h(k_h,1)    = mean(t_seg);
        freq_h_window(k_h,1) = bpm_h;

        for ii = idx_win
            dist = abs(ii - center_theoretical);
            w    = (halfN + 1) - dist;
            if w < 0, w = 0; end

            freq_weighted_sum_h(ii) = freq_weighted_sum_h(ii) + bpm_h * w;
            weight_sum_h(ii)        = weight_sum_h(ii)        + w;
        end
    end
end

%% ====== 1) FREQ_H POR AMOSTRA ======
freq_h_raw = nan(nt,1);
valid_h = weight_sum_h > 0;
freq_h_raw(valid_h) = freq_weighted_sum_h(valid_h) ./ weight_sum_h(valid_h);

%% ====== 2) SUAVIZAÇÃO H POR SEGMENTO ======
freq_h_smooth = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_h_smooth(idx_seg) = movmean(freq_h_raw(idx_seg), win_smooth_raw, 'omitnan');
end

%% ====== 3) INTERPOLAÇÃO H PARA FECHAR BURACOS ======
idx_good_h = ~isnan(freq_h_smooth);
if any(idx_good_h)
    freq_h_smooth(~idx_good_h) = interp1(t(idx_good_h), freq_h_smooth(idx_good_h), ...
                                         t(~idx_good_h),'linear','extrap');
end

%% ====== 4) SUAVIZAÇÃO FINAL H ======
freq_h_final = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_h_final(idx_seg) = movmean(freq_h_smooth(idx_seg), win_smooth_final, 'omitnan');
end

%% ======================================================================
%% ========================= AGORA BREATH ===============================
%% ======================================================================

%% ====== FREX PELO ESPECTRO DO BREATH ======
Nfft_br  = 2*nt - 1;
f_fft_br = linspace(0, srate, Nfft_br);
B        = fft(breath, Nfft_br);
[~, imax_br] = max(abs(B));
frex_br = f_fft_br(imax_br);

%% ====== DEFINIR h_br = X / frex_br ======
h_br = X / frex_br;
h_br = max(h_br, 0.1);
h_br = min(h_br, 5.0);

%% ========================== WAVELET & CONVOLUÇÃO (BREATH) ========================== %%
gaus_br = exp((-4*log(2)*t_aux.^2)/h_br^2) .* exp(1i*2*pi*frex_br*t_aux);

nGaus_br  = length(gaus_br);
nConv_br  = nGaus_br + nt - 1;

breathX = fft(breath, nConv_br);
gausX_br  = fft(gaus_br,  nConv_br);
gausX_br  = gausX_br ./ max(gausX_br);

convX_br     = breathX .* gausX_br;
convs_all_br = ifft(convX_br);

startIdx_br = floor(nGaus_br/2) + 1;
endIdx_br   = startIdx_br + nt - 1;
convs_br    = convs_all_br(startIdx_br:endIdx_br);

%% ====== FFT POR JANELAS (BREATH, SEGMENTADO) ======
freq_weighted_sum_br = zeros(nt,1);
weight_sum_br        = zeros(nt,1);

t_window_br      = [];
freq_br_window   = [];
k_br = 0;

for s = 1:nSeg
    seg_start = seg_starts(s);
    seg_end   = seg_ends(s);

    center_list = seg_start : hop : seg_end;

    for center_theoretical = center_list
        idx_theo = center_theoretical + rel_idx;
        idx_win  = idx_theo(idx_theo >= seg_start & idx_theo <= seg_end);

        if numel(idx_win) < 4
            continue;
        end

        t_seg = t(idx_win);
        x_seg = real(convs_br(idx_win));

        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end

        fs_seg = 1/mean(dt_seg);

        x_seg = x_seg - mean(x_seg);

        Xw   = fft(x_seg);
        Lwin = numel(x_seg);
        half = floor(Lwin/2);

        f_ax = (0:half-1) * fs_seg / Lwin;
        magX = abs(Xw(1:half));

        mask_br = (f_ax >= f_low_br) & (f_ax <= f_high_br);
        if ~any(mask_br), continue; end

        [~, imax_win_br] = max(magX(mask_br));
        f_est_vec_br = f_ax(mask_br);
        f_est_br     = f_est_vec_br(imax_win_br);

        bpm_br = f_est_br * 60;  % respiração em "bpm" (na prática, brpm)

        k_br = k_br + 1;
        t_window_br(k_br,1)    = mean(t_seg);
        freq_br_window(k_br,1) = bpm_br;

        for ii = idx_win
            dist = abs(ii - center_theoretical);
            w    = (halfN + 1) - dist;
            if w < 0, w = 0; end

            freq_weighted_sum_br(ii) = freq_weighted_sum_br(ii) + bpm_br * w;
            weight_sum_br(ii)        = weight_sum_br(ii)        + w;
        end
    end
end

%% ====== 1) FREQ_BR POR AMOSTRA ======
freq_br_raw = nan(nt,1);
valid_br = weight_sum_br > 0;
freq_br_raw(valid_br) = freq_weighted_sum_br(valid_br) ./ weight_sum_br(valid_br);

%% ====== 2) SUAVIZAÇÃO BR POR SEGMENTO ======
freq_br_smooth = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_br_smooth(idx_seg) = movmean(freq_br_raw(idx_seg), win_smooth_raw, 'omitnan');
end

%% ====== 3) INTERPOLAÇÃO BR PARA FECHAR BURACOS ======
idx_good_br = ~isnan(freq_br_smooth);
if any(idx_good_br)
    freq_br_smooth(~idx_good_br) = interp1(t(idx_good_br), freq_br_smooth(idx_good_br), ...
                                           t(~idx_good_br),'linear','extrap');
end

%% ====== 4) SUAVIZAÇÃO FINAL BR ======
freq_br_final = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_br_final(idx_seg) = movmean(freq_br_smooth(idx_seg), win_smooth_final, 'omitnan');
end

%% ====== GRÁFICOS ======

% HEART
figure;
plot(t_window_h, freq_h_window, 'o'); hold on;
plot(t, freq_h_raw, '.');
plot(t, freq_h_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('HR (bpm)');
title(sprintf('HEART - HR final — h_h = %.3f (X=0.6)', h_h));
grid on;
legend('Janela','Raw','Final');

% BREATH
figure;
plot(t_window_br, freq_br_window, 'o'); hold on;
plot(t, freq_br_raw, '.');
plot(t, freq_br_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('RR (bpm eq.)');
title(sprintf('BREATH - RR final — h_{br} = %.3f (X=0.6)', h_br));
grid on;
legend('Janela','Raw','Final');

%% ====== SALVAR RESULTADO ======
freq_h       = freq_h_raw;
freq_h_cont  = freq_h_final;
mean_h_cont  = freq_h_final;

freq_br      = freq_br_raw;
freq_br_cont = freq_br_final;
mean_br_cont = freq_br_final;

save('freq_results.mat', 't', ...
    'freq_h_raw','freq_h', ...
    'freq_h_cont','mean_h_cont', ...
    'freq_br_raw','freq_br', ...
    'freq_br_cont','mean_br_cont', ...
    't_window_h','freq_h_window', ...
    't_window_br','freq_br_window', ...
    'seg_starts','seg_ends', ...
    'N','overlap','hop', ...
    'win_smooth_raw','win_smooth_final', ...
    'frex_h','h_h','frex_br','h_br','X');
