clc; clear; close all;

%% ===================== PARÂMETROS CONFIGURÁVEIS ===================== %%
CSV_PATH  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases.csv';
POLAR_HR_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH2.txt';  % <<<< ADICIONADO

COL_T_MS  = 1;
COL_HEART = 4;

gap_thr = 0.5;   % limiar de buraco (s)

% ------------------ CWT (compatível) ------------------
WAVELET        = 'morse';   % mais compatível com versões antigas
VOICES_PER_OCT = 48;       % 12/16/24/32

% ------------------ Ridge com continuidade em Hz -----------------------
% Score = Pnorm(:,t) - LAMBDA_HZ2 * (f - f_prev)^2
RIDGE_LAMBDA_HZ2 = 500;      % 10, 40, 100, 250...

% ------------------ Correção de harmônico (2x) -------------------------
HARM_RATIO = 0.8;          % se energia em f/2 >= ratio*energia em f => troca pra f/2

% ------------------ Suavização temporal (sem cutoffs em Hz) ------------
MEDIAN_SMOOTH_SEC = 0;      % mediana móvel (mata spikes)
MEAN_SMOOTH_SEC   = 1.0;  % média móvel (alisa a curva)

FILL_NAN_WITHIN_SEG = true;

% ------------------ Comparação com Polar -------------------------------
DO_POLAR_COMPARE = true;
DT_GRID          = 0.1;     % igual aos outros
IGNORE_START_SEC = 0;       % igual aos outros

% ------------------ Plot ----------------------------------------------
PLOT_HR = true;

%% ========================== LEITURA ========================== %%
if ~isfile(CSV_PATH)
    error('CSV não encontrado: %s', CSV_PATH);
end

data  = readmatrix(CSV_PATH);
if size(data,2) < max(COL_T_MS, COL_HEART)
    error('CSV tem %d colunas; COL_T_MS=%d, COL_HEART=%d.', size(data,2), COL_T_MS, COL_HEART);
end

tpuro = data(:,COL_T_MS);
heart = data(:,COL_HEART);

t  = (tpuro - tpuro(1)) / 1000;
nt = numel(t);

fprintf('CSV OK: %d linhas | t(1)=%.0f ms | heart(1)=%.6f\n', size(data,1), tpuro(1), heart(1));

%% ==================== SEGMENTAÇÃO POR BURACOS ==================== %%
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

%% ====== dt e srate NOMINAIS (referência global) ======
dt_good = dt_all(dt_all > 0 & dt_all <= gap_thr);
if isempty(dt_good)
    dt = mean(dt_all(dt_all > 0));
else
    dt = mean(dt_good);
end
srate_nom = 1/dt;

%% ======================================================================
%% ===================== HR via CWT ridge (POR SEGMENTO) =================
%% ======================================================================

hr_hz_raw  = nan(nt,1);
hr_bpm_raw = nan(nt,1);

for s = 1:nSeg
    idx_seg = seg_starts(s):seg_ends(s);
    L = numel(idx_seg);
    if L < 64, continue; end

    t_seg = t(idx_seg);
    x_seg = heart(idx_seg);

    dt_seg = diff(t_seg);
    if any(dt_seg <= 0), continue; end
    fs = 1/mean(dt_seg);

    % remove DC (sem filtros)
    x_seg = x_seg - mean(x_seg);

    % --------- CWT compatível ---------
    [wt,f] = cwt(x_seg, fs, WAVELET, 'VoicesPerOctave',VOICES_PER_OCT);
    f = f(:);
    [f, ord] = sort(f,'ascend');
    wt = wt(ord,:);

    P = abs(wt).^2;                 % nF x L
    if isempty(P) || size(P,2) ~= L, continue; end

    % evita DC: zera apenas o 1º bin (sem cutoff em Hz)
    P(1,:) = 0;

    % --------- NORMALIZA por coluna ---------
    Pnorm = P ./ (max(P,[],1) + eps);

    % --------- ridge com continuidade em HZ ------------
    ridge_f = ridge_track_continuity_hz(Pnorm, f, RIDGE_LAMBDA_HZ2);

    % --------- correção de harmônico 2x ----------------
    ridge_f = harmonic_fix_half(P, f, ridge_f, HARM_RATIO);

    hr_hz  = ridge_f(:);
    hr_bpm = 60 * hr_hz;

    % salva bruto
    hr_hz_raw(idx_seg)  = hr_hz;
    hr_bpm_raw(idx_seg) = hr_bpm;

    % --------- suavização dentro do segmento -----------
    wMed = max(3, round(MEDIAN_SMOOTH_SEC * fs));
    wMean= max(3, round(MEAN_SMOOTH_SEC   * fs));

    tmp = hr_bpm_raw(idx_seg);

    tmp = movmedian(tmp, wMed, 'omitnan');

    if FILL_NAN_WITHIN_SEG
        tmp = fill_nan_by_interp(t_seg, tmp);
    end

    tmp = movmean(tmp, wMean, 'omitnan');

    hr_bpm_raw(idx_seg) = tmp;
end

%% ===================== FINAL =====================
hr_bpm_final = hr_bpm_raw;  % já suavizado por segmento

%% ===================== COMPARAÇÃO COM POLAR (IGUAL AOS OUTROS) =====================
if DO_POLAR_COMPARE
    if ~isfile(POLAR_HR_PATH)
        warning('Polar não encontrado: %s (pulando comparação)', POLAR_HR_PATH);
    else
        [t_polar, HR_polar] = read_txt_polar(POLAR_HR_PATH);
        
        % remove duplicatas
        [t_polar, ia] = unique(t_polar, 'stable');
        HR_polar = HR_polar(ia);

        HR_polar_mean = mean(HR_polar, 'omitnan');
        fprintf('HR Polar médio: %.4f bpm | f média: %.6f Hz\n', HR_polar_mean, HR_polar_mean/60);

        t_radar  = t;
        HR_radar = hr_bpm_final;

        okR = isfinite(t_radar) & isfinite(HR_radar);
        okP = isfinite(t_polar) & isfinite(HR_polar);

        if nnz(okR) < 2 || nnz(okP) < 2
            warning('Poucos pontos válidos Radar x Polar (radar=%d, polar=%d).', nnz(okR), nnz(okP));
        else
            t_start = max(t_radar(find(okR,1,'first')), t_polar(find(okP,1,'first')));
            t_end   = min(t_radar(find(okR,1,'last')),  t_polar(find(okP,1,'last')));

            if ~(isfinite(t_start) && isfinite(t_end)) || (t_end <= t_start)
                warning('Sem interseção de tempo válida entre Radar e Polar.');
            else
                t_common = (t_start : DT_GRID : t_end)';

                HR_radar_grid = interp1(t_radar(okR), HR_radar(okR), t_common, 'linear');
                HR_polar_grid = interp1(t_polar(okP), HR_polar(okP), t_common, 'linear');

                mask_time = (t_common >= (t_start + IGNORE_START_SEC));

                HR_ref  = HR_polar_grid(mask_time);
                HR_test = HR_radar_grid(mask_time);

                err = HR_test - HR_ref;

                MAE  = mean(abs(err), 'omitnan');
                RMSE = sqrt(mean(err.^2, 'omitnan'));
                R    = corr(HR_test, HR_ref, 'rows', 'complete');

                fprintf('\n=== Radar x Polar (dt_grid=%.2fs | ignore=%.1fs) ===\n', DT_GRID, IGNORE_START_SEC);
                fprintf('MAE : %.4f bpm\n', MAE);
                fprintf('RMSE: %.4f bpm\n', RMSE);
                fprintf('Corr: %.4f\n', R);
            end
        end
    end
end

%% ===================== PLOT FINAL =====================
if PLOT_HR
    figure('Color','w');
    plot(t, hr_bpm_final, 'k', 'LineWidth', 1.4); hold on;
    grid on;
    xlabel('Tempo (s)');
    ylabel('HR (bpm)');
    title(sprintf('HEART | HR via CWT ridge (contínuo) | srate_nom=%.2f Hz | nSeg=%d', srate_nom, nSeg));

    % overlay Polar (se existir)
    if DO_POLAR_COMPARE && isfile(POLAR_HR_PATH)
        [t_polar, HR_polar] = read_txt_polar(POLAR_HR_PATH);
        [t_polar, ia] = unique(t_polar, 'stable');
        HR_polar = HR_polar(ia);
        plot(t_polar, HR_polar, 'm.-', 'LineWidth', 1, 'MarkerSize', 8);
        legend('Radar (final)', 'Polar', 'Location','best');
    else
        legend('Radar (final)', 'Location','best');
    end
end

%% ===================== (SEM SAVE) =====================
% Você pediu "sem save nenhum", então removi o save().

%% ===================== FUNÇÕES LOCAIS =====================
function [cfs,f] = get_cwt_compat(x, fs, wname, vpo)
    L = numel(x);

    if exist('cwtfilterbank','file')
        try
            fb = cwtfilterbank('SignalLength', L, ...
                               'SamplingFrequency', fs, ...
                               'VoicesPerOctave', vpo, ...
                               'Wavelet', wname);
        catch
            fb = cwtfilterbank('SignalLength', L, ...
                               'SamplingFrequency', fs, ...
                               'VoicesPerOctave', vpo, ...
                               'Wavelet', 'morse');
        end

        try
            [cfs,f] = cwt(fb, x); return;
        catch
        end
        try
            [cfs,f] = cwt(x, 'FilterBank', fb); return;
        catch
        end
        try
            [cfs,f] = fb.cwt(x); return;
        catch
        end

        error('cwtfilterbank existe, mas não encontrei forma compatível de aplicar a CWT.');
    end

    try
        [cfs,f] = cwt(x, fs); return;
    catch
    end
    try
        [cfs,f] = cwt(x, 1/fs); return;
    catch
    end

    error('Sem suporte a CWT nesta instalação.');
end

function ridge_f = ridge_track_continuity_hz(Pnorm, f, lambda_hz2)
    [nF,L] = size(Pnorm);
    ridge_f = nan(1,L);

    [~, i0] = max(Pnorm(:,1));
    ridge_f(1) = f(i0);
    fprev = ridge_f(1);

    for k = 2:L
        penal = lambda_hz2 * (f - fprev).^2;
        score = Pnorm(:,k) - penal;
        [~, ik] = max(score);
        ridge_f(k) = f(ik);
        fprev = ridge_f(k);
    end
end

function ridge_f = harmonic_fix_half(P, f, ridge_f, ratio)
    L = size(P,2);
    for k = 1:L
        fk = ridge_f(k);
        if isnan(fk) || fk <= 0, continue; end
        [~, i1] = min(abs(f - fk));
        [~, i2] = min(abs(f - fk/2));
        p1 = P(i1,k);
        p2 = P(i2,k);
        if p2 >= ratio * p1
            ridge_f(k) = fk/2;
        end
    end
end

function x = fill_nan_by_interp(t, x)
    ok = ~isnan(x);
    if nnz(ok) >= 2
        x(~ok) = interp1(t(ok), x(ok), t(~ok), 'linear', 'extrap');
    end
end

function [t_sec, HR] = read_txt_polar(p)
    fid = fopen(p,'r');
    if fid == -1, error('Não foi possível abrir o arquivo Polar: %s', p); end
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
