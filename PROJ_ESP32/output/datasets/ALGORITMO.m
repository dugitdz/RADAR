clc; clear; close all;

%% ===================== PARÂMETROS CONFIGURÁVEIS ===================== %%
N        = 1024;     % janela (amostras) p/ FFT em janelas
overlap  = 0.8;
hop      = max(1, round(N*(1-overlap)));

win_smooth_raw   = 6400;
win_smooth_final = 6400;

gap_thr = 0.5;       % limiar de buraco (s)

% Wavelet: h = X/frex
X = 1.2;

% ===== FREX (estilo antigo) =====
USE_FREX_EST = false;        % true: estima frex do heart; false: usa frex_fix
frex_fix     = 1.4;          % Hz (90 bpm) fallback
HR_BAND_HZ   = [];      % faixa plausível p/ HR (Hz). Use [] pra não restringir

%% ========================== SAÍDA ========================== %%
OUT_DIR = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\datasets\Mendeley\output';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ========================== LEITURA (CSV OU MAT) ========================== %%
INPUT_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\datasets\Mendeley\subject8_phaseHR.mat';
[tpuro, heart, baseName] = load_phase_for_algorithm(INPUT_PATH);

%% ====== TEMPO ======
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

%% ====== dt e srate NOMINAIS ======
dt_good = dt_all(dt_all > 0 & dt_all <= gap_thr);
if isempty(dt_good)
    dt = mean(dt_all(dt_all > 0));
else
    dt = mean(dt_good);
end
srate = 1/dt;

%% ======================================================================
%% =============== FREX (ESTILO ANTIGO: Nfft=2*L-1) ======================
%% ======================================================================
if USE_FREX_EST
    seg_frex = nan(nSeg,1);
    seg_w    = nan(nSeg,1);

    for s = 1:nSeg
        idx_seg = seg_starts(s):seg_ends(s);
        Ls = numel(idx_seg);
        if Ls < 32, continue; end   % muito pequeno => lixo

        t_seg = t(idx_seg);
        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end
        fs_seg = 1/mean(dt_seg);

        xs = heart(idx_seg);
        xs = xs - mean(xs);

        Nfft = 2*Ls - 1;
        H    = fft(xs, Nfft);

        half = floor(Nfft/2) + 1;
        magH = abs(H(1:half));
        f_ax = (0:half-1) * (fs_seg / Nfft);

        magH(1) = 0; % remove DC

        if ~isempty(HR_BAND_HZ)
            band = (f_ax >= HR_BAND_HZ(1)) & (f_ax <= HR_BAND_HZ(2));
            if ~any(band), continue; end
            mag2 = magH; mag2(~band) = 0;
            [~, imax] = max(mag2);
        else
            [~, imax] = max(magH);
        end

        frex = f_ax(imax);
        if ~isfinite(frex) || frex <= 0, continue; end

        seg_frex(s) = frex;
        seg_w(s)    = Ls;
    end

    ok = isfinite(seg_frex) & isfinite(seg_w) & seg_w > 0;
    if any(ok)
        frex_h = sum(seg_frex(ok).*seg_w(ok)) / sum(seg_w(ok));
    else
        frex_h = frex_fix;
    end
else
    frex_h = frex_fix;
end

fprintf('frex_h = %.6f Hz (%.2f bpm) | srate_nom=%.3f Hz\n', frex_h, 60*frex_h, srate);

%% ======================================================================
%% ===================== HEART (WAVELET + JANELAS) =======================
%% ======================================================================

%% ====== h_h = X / frex_h ======
h_h = X / frex_h;
h_h = max(h_h, 0.1);
h_h = min(h_h, 5.0);

%% ====== WAVELET & CONVOLUÇÃO (POR SEGMENTO) ======
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
    gausX  = gausX ./ (max(abs(gausX)) + eps);

    convX  = heartX .* gausX;
    y_all  = ifft(convX);

    startIdx = floor(nG/2) + 1;
    endIdx   = startIdx + L - 1;

    convs_h(idx_seg) = y_all(startIdx:endIdx);
end

%% ====== FFT POR JANELAS (HEART) ======
halfN   = floor(N/2);
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

        % só aceita janela COMPLETA
        if any(idx_theo < seg_start) || any(idx_theo > seg_end)
            continue;
        end
        idx_win = idx_theo;

        t_seg = t(idx_win);
        x_seg = real(convs_h(idx_win));
        if any(isnan(x_seg)), continue; end

        % fs por janela
        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end
        fs_seg = 1/mean(dt_seg);

        x_seg = x_seg - mean(x_seg);

        Lwin = numel(x_seg);
        if Lwin < 8, continue; end

        % >>> FFT "antigo": Nfft = 2*Lwin - 1
        Nfft = 2*Lwin - 1;
        Xw   = fft(x_seg, Nfft);

        half = floor(Nfft/2) + 1;
        f_ax = (0:half-1) * fs_seg / Nfft;
        magX = abs(Xw(1:half));

        magX(1) = 0; % remove DC

        % opcional: restringir busca em faixa plausível (evita respiração)
        if ~isempty(HR_BAND_HZ)
            band = (f_ax >= HR_BAND_HZ(1)) & (f_ax <= HR_BAND_HZ(2));
            if ~any(band), continue; end
            mag2 = magX; mag2(~band) = 0;
            [~, imax_win] = max(mag2);
        else
            [~, imax_win] = max(magX);
        end

        % ===================== CORREÇÃO DO "500 bpm" =====================
        % refinamento parabólico blindado:
        % 1) só se for pico real (b>a e b>c)
        % 2) denom relativo (evita explode quando espectro é plano)
        % 3) limita |p| <= 0.5 (parábola só faz sentido entre bins adjacentes)
        p = 0;
        if imax_win > 1 && imax_win < numel(magX)
            a = magX(imax_win-1);
            b = magX(imax_win);
            c = magX(imax_win+1);

            if (b > a) && (b > c)
                denom = (a - 2*b + c);
                scale = max([a b c 1]);          % evita scale=0
                if abs(denom) > 1e-6 * scale     % << blindagem forte
                    p = 0.5*(a - c)/denom;
                    p = max(min(p, 0.5), -0.5);  % << limite crítico
                end
            end
        end

        f_est_h = ((imax_win-1) + p) * fs_seg / Nfft;

        % garante que NÃO sai da banda por causa do p/ruído
        if ~isempty(HR_BAND_HZ)
            f_est_h = min(max(f_est_h, HR_BAND_HZ(1)), HR_BAND_HZ(2));
        end

        bpm_h   = 60 * f_est_h;
        % =================================================================

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

%% ====== GRÁFICO ======
figure;
plot(t_window_h, freq_h_window, 'o'); hold on;
plot(t, freq_h_raw, '.');
plot(t, freq_h_final, 'k', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('HR (bpm)');
title(sprintf('HEART | frex_h=%.3f Hz (%.1f bpm) | h_h=%.3f | srate=%.2f Hz | N=%d', ...
    frex_h, frex_h*60, h_h, srate, N));
grid on;
legend('Janela','Raw','Final');

%% ====== SALVAR ======
freq_h       = freq_h_raw;
freq_h_cont  = freq_h_final;
mean_h_cont  = freq_h_final;

out_file = fullfile(OUT_DIR, [baseName '_freq_results.mat']);

save(out_file, 't', ...
    'freq_h_raw','freq_h','freq_h_cont','mean_h_cont', ...
    't_window_h','freq_h_window', ...
    'seg_starts','seg_ends', ...
    'N','overlap','hop', ...
    'win_smooth_raw','win_smooth_final', ...
    'gap_thr','dt','srate', ...
    'frex_h','h_h','X', ...
    'USE_FREX_EST','frex_fix','HR_BAND_HZ');

fprintf('\nSalvo em: %s\n', out_file);

%% ========================== FUNÇÕES LOCAIS ==========================
function [tpuro_ms, heart, baseName] = load_phase_for_algorithm(path_in)
    [~, baseName, ext] = fileparts(path_in);
    ext = lower(ext);

    if strcmp(ext,'.csv')
        data = readmatrix(path_in);
        tpuro_ms = data(:,1);
        phase    = data(:,3);
        heart = wrap_phase(phase);
        return;
    end

    if strcmp(ext,'.mat')
        S = load(path_in);
        if ~isfield(S,'DS'), error('MAT precisa conter struct DS.'); end
        DS = S.DS;

        if ~isfield(DS,'t'), error('DS precisa ter DS.t'); end
        tpuro_ms = DS.t(:) * 1000;

        % blindagem: recalcula fase do radar se tiver I/Q
        if isfield(DS,'radar24I') && isfield(DS,'radar24Q')
            I = DS.radar24I(:);
            Q = DS.radar24Q(:);
            phase = unwrap(atan2(Q, I));
        elseif isfield(DS,'phase')
            phase = DS.phase(:);
        else
            error('DS precisa ter radar24I/radar24Q ou DS.phase.');
        end

        heart = wrap_phase(phase);
        return;
    end

    error('Entrada inválida: use .csv ou .mat (DS).');
end

function ph_wrapped = wrap_phase(ph)
    ph = ph(:);

    ok = isfinite(ph);
    if ~all(ok)
        ph(~ok) = interp1(find(ok), ph(ok), find(~ok), 'linear', 'extrap');
    end

    ph_wrapped = mod(ph + pi, 2*pi) - pi; % [-pi,pi]
    ph_wrapped = ph_wrapped - mean(ph_wrapped);
end
