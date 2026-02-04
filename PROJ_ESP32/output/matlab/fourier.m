clc; clear; close all;

%% ===================== PARÂMETROS CONFIGURÁVEIS ===================== %%
N        = 16;
overlap  = 0.8;
hop      = max(1, round(N*(1-overlap)));

win_smooth_raw   = 128;
win_smooth_final = 128;

gap_thr = 0.5;   % limiar de buraco (s)

% Wavelet: h = X/frex
X = 0.6;         % fixo como no seu código

% ===== PLOT espectro por segmento (cada um em uma figura) =====
PLOT_SEG_SPECTRA = true;
XLIM_HZ = [0 10];    % ajuste se quiser (ex.: [0 5], [] pra não limitar)

%% ========================== LEITURA ========================== %%
data   = readmatrix('C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phase.csv');
tpuro  = data(:,1);
breath = data(:,3);   %#ok<NASGU>  % (desativado mais abaixo)
heart  = data(:,4);

%% ====== TEMPO (sem srate ainda) ======
t  = (tpuro - tpuro(1)) / 1000;
nt = numel(t);

%% ====== SEGMENTAÇÃO POR BURACOS ======
dt_all   = diff(t);
break_ix = find(dt_all > gap_thr);

if isempty(break_ix)
    seg_starts = 1;
    seg_ends   = nt;
else
    seg_starts = [1; break_ix + 1];
    seg_ends   = [break_ix; nt];
end
nSeg = numel(seg_starts);

%% ====== dt e srate NOMINAIS (só com dt "bons", sem mediana) ======
dt_good = dt_all(dt_all > 0 & dt_all <= gap_thr);

if isempty(dt_good)
    dt = mean(dt_all(dt_all > 0));   % fallback extremo
else
    dt = mean(dt_good);              % <<< sem mediana
end
srate = 1/dt;                        % só referência global

%% ======================================================================
%% ============ FREX "IGUAL AO ANTIGO", MAS POR SEGMENTO (HEART) =========
%% ============ + PLOT DO ESPECTRO (1 FIGURA POR SEGMENTO) ===============
%% ======================================================================

fprintf('\n=== HEART: pico por segmento (FFT estilo antigo) ===\n');

seg_fpk = nan(nSeg,1);

wsum = 0;
fsum = 0;

for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    L = numel(idx_seg);
    if L < N, continue; end

    t_seg = t(idx_seg);
    x_seg = heart(idx_seg);

    dt_seg = diff(t_seg);
    if any(dt_seg <= 0), continue; end
    fs_seg = 1/mean(dt_seg);

    x_seg = x_seg - mean(x_seg);  % reduz DC

    % === "antigo", só que com L e fs_seg ===
    Nfft  = 2*L - 1;
    f_fft = linspace(0, fs_seg, Nfft);
    H     = fft(x_seg, Nfft);

    halfSpec = floor(Nfft/2);
    if halfSpec < 3, continue; end

    magH = abs(H(1:halfSpec));

    % pico (ignorando DC)
    [~, imax] = max(magH(2:end));
    imax = imax + 1;
    fpk = f_fft(imax);
    seg_fpk(s) = fpk;

    % média ponderada básica (peso = L)
    w = L;
    fsum = fsum + fpk * w;
    wsum = wsum + w;

    fprintf('seg %3d | L=%5d | fs=%.3f Hz | fpk=%.5f Hz (%.1f bpm)\n', s, L, fs_seg, fpk, 60*fpk);

    % ===== PLOT (cada segmento em uma figura) =====
    if PLOT_SEG_SPECTRA
        magPlot = magH;
        if max(magPlot) > 0
            magPlot = magPlot ./ max(magPlot);  % normaliza por segmento
        end

        figure;
        plot(f_fft(1:halfSpec), magPlot, '.', 'MarkerSize', 6); % <<< SÓ PONTOS
        grid on;
        xlabel('Frequência (Hz)');
        ylabel('Magnitude (normalizada)');
        title(sprintf('HEART | Segmento %d | L=%d | fs=%.3f Hz | fpk=%.5f Hz (%.1f bpm)', ...
            s, L, fs_seg, fpk, 60*fpk));

        % marca o pico
        xline(fpk, '--');

        if ~isempty(XLIM_HZ)
            xlim(XLIM_HZ);
        end
    end
end

if wsum > 0
    frex_h = fsum/wsum;
else
    error('Não foi possível estimar frex_h: nenhum segmento válido.');
end

fprintf('-----------------------------------------------\n');
fprintf('frex_h (média ponderada) = %.6f Hz (%.2f bpm)\n\n', frex_h, 60*frex_h);

%% ======================================================================
%% ===================== HEART (WAVELET + JANELAS) =======================
%% ======================================================================

%% ====== h_h = X / frex_h ======
h_h = X / frex_h;
h_h = max(h_h, 0.1);
h_h = min(h_h, 5.0);

%% ====== WAVELET & CONVOLUÇÃO (POR SEGMENTO) - HEART ======
convs_h = nan(nt,1);

for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    L = numel(idx_seg);
    if L < 8, continue; end

    t_seg = t(idx_seg);
    dt_seg = diff(t_seg);
    if any(dt_seg <= 0), continue; end
    dt_loc = mean(dt_seg);

    idx_loc = (1:L) - ceil(L/2);
    t_aux   = idx_loc(:) * dt_loc;

    gaus_h = exp((-4*log(2)*t_aux.^2)/h_h^2) .* exp(1i*2*pi*frex_h*t_aux);

    nG = length(gaus_h);
    nConv = nG + L - 1;

    x = heart(idx_seg);
    heartX = fft(x, nConv);
    gausX  = fft(gaus_h, nConv);
    gausX  = gausX ./ max(gausX);

    convX  = heartX .* gausX;
    y_all  = ifft(convX);

    startIdx = floor(nG/2) + 1;
    endIdx   = startIdx + L - 1;

    convs_h(idx_seg) = y_all(startIdx:endIdx);
end

%% ====== FFT POR JANELAS (HEART) ======
halfN   = N/2;
rel_idx = -halfN : (halfN-1);

freq_weighted_sum_h = zeros(nt,1);
weight_sum_h        = zeros(nt,1);

t_window_h    = [];
freq_h_window = [];
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
        if any(isnan(x_seg)), continue; end

        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end
        fs_seg = 1/mean(dt_seg);

        x_seg = x_seg - mean(x_seg);

        Xw   = fft(x_seg);
        Lwin = numel(x_seg);
        half = floor(Lwin/2);
        if half < 3, continue; end

        f_ax = (0:half-1) * fs_seg / Lwin;
        magX = abs(Xw(1:half));

        [~, imax_win] = max(magX(2:end));  % ignora DC
        imax_win = imax_win + 1;
        f_est_h  = f_ax(imax_win);

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

freq_h_raw = nan(nt,1);
valid_h = weight_sum_h > 0;
freq_h_raw(valid_h) = freq_weighted_sum_h(valid_h) ./ weight_sum_h(valid_h);

freq_h_smooth = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_h_smooth(idx_seg) = movmean(freq_h_raw(idx_seg), win_smooth_raw, 'omitnan');
end

idx_good_h = ~isnan(freq_h_smooth);
if any(idx_good_h)
    freq_h_smooth(~idx_good_h) = interp1(t(idx_good_h), freq_h_smooth(idx_good_h), ...
                                         t(~idx_good_h),'linear','extrap');
end

freq_h_final = nan(nt,1);
for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    freq_h_final(idx_seg) = movmean(freq_h_smooth(idx_seg), win_smooth_final, 'omitnan');
end

%% ====== GRÁFICO HEART ======
figure;
plot(t_window_h, freq_h_window, 'o'); hold on;
plot(t, freq_h_raw, '.');
plot(t, freq_h_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('HR (bpm)');
title(sprintf('HEART | frex_h=%.3f Hz (%.1f bpm) | h_h=%.3f | srate_nom=%.2f Hz', ...
    frex_h, frex_h*60, h_h, srate));
grid on;
legend('Janela','Raw','Final');

%% ======================================================================
%% ===================== BREATH (DESATIVADO) =============================
%% ======================================================================
if false
    % (breath desativado)
end

%% ====== SALVAR RESULTADO ======
freq_h       = freq_h_raw;
freq_h_cont  = freq_h_final;
mean_h_cont  = freq_h_final;

save('freq_results.mat', 't', ...
    'freq_h_raw','freq_h','freq_h_cont','mean_h_cont', ...
    't_window_h','freq_h_window', ...
    'seg_starts','seg_ends', ...
    'N','overlap','hop', ...
    'win_smooth_raw','win_smooth_final', ...
    'gap_thr','dt','srate', ...
    'frex_h','h_h','X');
