clc; clear; close all;

%% ===================== PARÂMETROS CONFIGURÁVEIS ===================== %%
N        = 128;
overlap  = 0.8;
hop      = max(1, round(N*(1-overlap)));

win_smooth_raw   = 128;
win_smooth_final = 128;

% ====== CORTES DO PASSA-FAIXA (AJUSTE AQUI) ======
% HEART em bpm:
heart_bp_bpm = [50 150];   % <<<<--- pedido: 60;150

% BREATH em brpm (sugestão default):
breath_bp_brpm = [6 30];   % ajuste como quiser

% filtro
bp_order_h  = 2;           % ordem Butterworth (heart)
bp_order_br = 4;           % ordem Butterworth (breath)
min_f_hz    = 0.05;        % evita encostar em 0 Hz
guard_hz    = 0.05;        % margem p/ não encostar em Nyquist

%% ===================== PATHS ===================== %%
PHASES_PATH    = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases.csv';
polar_HR_path  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH2.txt';

%% ========================== LEITURA ========================== %%
data   = readmatrix(PHASES_PATH);
tpuro  = data(:,1);
breath = data(:,3);
heart  = data(:,4);

%% ====== TEMPO ======
t  = (tpuro - tpuro(1)) / 1000;
dt = mean(diff(t));
srate = 1/dt; %#ok<NASGU>
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

%% ===================== FILTRO PASSA-FAIXA FIXO (HEART) POR SEGMENTO ===================== %%
convs_h = nan(nt,1);

% cortes em Hz
heart_bp_hz = heart_bp_bpm / 60;

for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    t_seg   = t(idx_seg);
    x_seg   = heart(idx_seg);

    if numel(idx_seg) < 8
        continue;
    end

    dt_seg = diff(t_seg);
    if any(dt_seg <= 0), continue; end
    fs_seg = 1/mean(dt_seg);

    f1 = max(min_f_hz, heart_bp_hz(1));
    f2 = min((fs_seg/2) - guard_hz, heart_bp_hz(2));

    if ~(isfinite(f1) && isfinite(f2) && f2 > f1)
        continue;
    end

    x_seg = x_seg - mean(x_seg);

    Wn = [f1 f2] / (fs_seg/2);
    if any(Wn <= 0) || any(Wn >= 1)
        continue;
    end

    [b,a] = butter(bp_order_h, Wn, 'bandpass');

    try
        y_seg = filtfilt(b,a,x_seg);
    catch
        continue;
    end

    convs_h(idx_seg) = y_seg;
end

%% ===================== FILTRO PASSA-FAIXA FIXO (BREATH) POR SEGMENTO ===================== %%
convs_br = nan(nt,1);

breath_bp_hz = breath_bp_brpm / 60;

for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    t_seg   = t(idx_seg);
    x_seg   = breath(idx_seg);

    if numel(idx_seg) < 8
        continue;
    end

    dt_seg = diff(t_seg);
    if any(dt_seg <= 0), continue; end
    fs_seg = 1/mean(dt_seg);

    f1 = max(min_f_hz, breath_bp_hz(1));
    f2 = min((fs_seg/2) - guard_hz, breath_bp_hz(2));

    if ~(isfinite(f1) && isfinite(f2) && f2 > f1)
        continue;
    end

    x_seg = x_seg - mean(x_seg);

    Wn = [f1 f2] / (fs_seg/2);
    if any(Wn <= 0) || any(Wn >= 1)
        continue;
    end

    [b,a] = butter(bp_order_br, Wn, 'bandpass');

    try
        y_seg = filtfilt(b,a,x_seg);
    catch
        continue;
    end

    convs_br(idx_seg) = y_seg;
end

%% ===================== FFT POR JANELAS (HEART, SEGMENTADO) ===================== %%
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
        x_seg = convs_h(idx_win);

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

        % pico dominante (ignora DC)
        [~, imax_win] = max(magX(2:end));
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

%% ===================== FFT POR JANELAS (BREATH, SEGMENTADO) ===================== %%
freq_weighted_sum_br = zeros(nt,1);
weight_sum_br        = zeros(nt,1);

t_window_br    = [];
freq_br_window = [];
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
        x_seg = convs_br(idx_win);

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

        [~, imax_win_br] = max(magX(2:end));
        imax_win_br = imax_win_br + 1;
        f_est_br    = f_ax(imax_win_br);

        bpm_br = f_est_br * 60;

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

%% ===================== POLAR: LEITURA (DO SEU JEITO) + ERROS ===================== %%
[t_polar, HR_polar] = read_txt_polar(polar_HR_path);

[t_polar, ia] = unique(t_polar, 'stable');
HR_polar = HR_polar(ia);

HR_polar_mean = mean(HR_polar(isfinite(HR_polar)));
f_polar_mean  = HR_polar_mean / 60; %#ok<NASGU>

% alinha (t_polar começa em 0) com eixo do radar (t também começa em 0)
t_polar_aligned = t_polar + t(1);

% recorta ao intervalo do radar
in = isfinite(t_polar_aligned) & isfinite(HR_polar) & ...
     (t_polar_aligned >= t(1)) & (t_polar_aligned <= t(end));

t_polar_aligned = t_polar_aligned(in);
HR_polar_in     = HR_polar(in);

% interpola radar nos tempos da Polar (AQUI NASCEM AS VARIÁVEIS)
HR_raw_at_polar   = interp1(t, freq_h_raw,   t_polar_aligned, 'linear', 'extrap');
HR_final_at_polar = interp1(t, freq_h_final, t_polar_aligned, 'linear', 'extrap');

% erros
err_raw   = HR_raw_at_polar   - HR_polar_in;
err_final = HR_final_at_polar - HR_polar_in;

% métricas compatíveis (sem omitnan)
v_raw = isfinite(err_raw)   & isfinite(HR_raw_at_polar)   & isfinite(HR_polar_in);
v_fin = isfinite(err_final) & isfinite(HR_final_at_polar) & isfinite(HR_polar_in);

metrics = struct();

if any(v_raw)
    e  = err_raw(v_raw);
    gt = HR_polar_in(v_raw);
    pr = HR_raw_at_polar(v_raw);

    metrics.raw.MAE  = mean(abs(e));
    metrics.raw.RMSE = sqrt(mean(e.^2));
    metrics.raw.BIAS = mean(e);
    metrics.raw.MAPE = mean(abs(e)./max(gt,1e-6))*100;
    metrics.raw.CORR = corr(pr(:), gt(:));
else
    metrics.raw.MAE  = NaN; metrics.raw.RMSE = NaN; metrics.raw.BIAS = NaN;
    metrics.raw.MAPE = NaN; metrics.raw.CORR = NaN;
end

if any(v_fin)
    e  = err_final(v_fin);
    gt = HR_polar_in(v_fin);
    pr = HR_final_at_polar(v_fin);

    metrics.final.MAE  = mean(abs(e));
    metrics.final.RMSE = sqrt(mean(e.^2));
    metrics.final.BIAS = mean(e);
    metrics.final.MAPE = mean(abs(e)./max(gt,1e-6))*100;
    metrics.final.CORR = corr(pr(:), gt(:));
else
    metrics.final.MAE  = NaN; metrics.final.RMSE = NaN; metrics.final.BIAS = NaN;
    metrics.final.MAPE = NaN; metrics.final.CORR = NaN;
end

fprintf('\n=== ERROS vs POLAR (HR) ===\n');
fprintf('RAW  : MAE=%.3f bpm | RMSE=%.3f bpm | BIAS=%.3f bpm | MAPE=%.2f%% | CORR=%.3f\n', ...
    metrics.raw.MAE, metrics.raw.RMSE, metrics.raw.BIAS, metrics.raw.MAPE, metrics.raw.CORR);
fprintf('FINAL: MAE=%.3f bpm | RMSE=%.3f bpm | BIAS=%.3f bpm | MAPE=%.2f%% | CORR=%.3f\n\n', ...
    metrics.final.MAE, metrics.final.RMSE, metrics.final.BIAS, metrics.final.MAPE, metrics.final.CORR);

%% ===================== PLOTS POLAR x RADAR + ERRO ===================== %%
figure;
plot(t_polar_aligned, HR_polar_in, 'k', 'LineWidth', 1.5); hold on;
plot(t_polar_aligned, HR_raw_at_polar, '.');
plot(t_polar_aligned, HR_final_at_polar, 'LineWidth', 1.2);
xlabel('Tempo (s)'); ylabel('HR (bpm)');
title('Comparação HR: Polar vs Radar (raw/final)');
grid on; legend('Polar','Radar@Polar (raw)','Radar@Polar (final)');

figure;
plot(t_polar_aligned, err_raw, '.'); hold on;
plot(t_polar_aligned, err_final, 'LineWidth', 1.2);
yline(0,'--');
xlabel('Tempo (s)'); ylabel('Erro (bpm)');
title('Erro HR (Radar - Polar)');
grid on; legend('Erro raw','Erro final','0');

%% ====== GRÁFICOS (SEUS) ======
figure;
plot(t_window_h, freq_h_window, 'o'); hold on;
plot(t, freq_h_raw, '.');
plot(t, freq_h_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('HR (bpm)');
title(sprintf('HEART | Bandpass %g-%g bpm | ordem %d', heart_bp_bpm(1), heart_bp_bpm(2), bp_order_h));
grid on;
legend('Janela','Raw','Final');

figure;
plot(t_window_br, freq_br_window, 'o'); hold on;
plot(t, freq_br_raw, '.');
plot(t, freq_br_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('RR (brpm)');
title(sprintf('BREATH | Bandpass %g-%g brpm | ordem %d', breath_bp_brpm(1), breath_bp_brpm(2), bp_order_br));
grid on;
legend('Janela','Raw','Final');

%% ===================== FUNÇÃO LOCAL (POLAR) ===================== %%
function [t_sec, HR] = read_txt_polar(p)
    fid = fopen(p,'r');
    if fid == -1, error('Não foi possível abrir o arquivo Polar.'); end
    C = textscan(fid,'%s %f %*[^\n]','Delimiter',';','HeaderLines',1);
    fclose(fid);

    ts_str = C{1};
    HR     = C{2};

    try
        t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
    catch
        t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    end

    t_sec = seconds(t_dt - t_dt(1));
end
