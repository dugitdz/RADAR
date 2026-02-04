clc; clear; close all;

%% ===================== PATHS =====================
CSV_PATH       = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases.csv';
POLAR_HR_PATH  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH2.txt';

COL_T_MS  = 1;
COL_HEART = 2;

gap_thr = 0.5;   % buraco (s)

%% ===================== FIXOS (fora do grid) =====================
FIX.FILL_NAN_WITHIN_SEG = false;

% suavização (movmean) -> FIXO (fora do grid, como você pediu)
FIX.MEAN_SMOOTH_SEC = 128.0;

% comparação Polar -> FIXO (fora do grid, como você pediu)
FIX.DO_POLAR_COMPARE = true;
FIX.DT_GRID          = 0.1;
FIX.IGNORE_START_SEC = 0;

% plot final do melhor (opcional)
FIX.PLOT_BEST = true;

%% ===================== GRID (tudo que varia) =====================
% CWT
GRID.WAVELET        = {'bpum','amor','morse'};   % tenta nessa ordem
GRID.VOICES_PER_OCT = [16];

% filtro (NOVO)
GRID.BP_ORDER = [2 6];
GRID.BP_WI    = [0.5 0.7];
GRID.BP_WF    = [2.5 3.5];

% ridge contínuo (Hz)
GRID.RIDGE_LAMBDA_HZ2 = [250];

% harmônico
GRID.HARM_RATIO = [0.75 0.85 0.95];

% suavização robusta (mediana)
GRID.MEDIAN_SMOOTH_SEC = [0 2 5];

PRINT_EVERY = 20;

%% ===================== CHECAGENS =====================
if ~isfile(CSV_PATH),  error('CSV não encontrado: %s', CSV_PATH); end
if ~isfile(POLAR_HR_PATH)
    warning('Polar não encontrado: %s (métricas ficarão NaN)', POLAR_HR_PATH);
end

%% ===================== ENUMERA GRID =====================
[A,B,C,D,E,F,G] = ndgrid( ...
    1:numel(GRID.WAVELET), ...
    GRID.VOICES_PER_OCT, ...
    GRID.BP_ORDER, ...
    GRID.BP_WI, ...
    GRID.BP_WF, ...
    GRID.RIDGE_LAMBDA_HZ2, ...
    GRID.HARM_RATIO);

% MEDIAN_SMOOTH_SEC entra como outro eixo (pra não explodir ndgrid acima)
MEDV = GRID.MEDIAN_SMOOTH_SEC(:);
comb0 = [A(:) B(:) C(:) D(:) E(:) F(:) G(:)];
n0 = size(comb0,1);
nComb = n0 * numel(MEDV);

fprintf('Total de combinações do grid: %d\n', nComb);

%% ===================== PREALOC RESULTADOS =====================
res = struct();
res.wavelet_id   = zeros(nComb,1);
res.vpo          = zeros(nComb,1);
res.bp_order     = zeros(nComb,1);
res.bp_wi        = zeros(nComb,1);
res.bp_wf        = zeros(nComb,1);
res.lambda_hz2   = zeros(nComb,1);
res.harm_ratio   = zeros(nComb,1);
res.med_sec      = zeros(nComb,1);

res.MAE  = nan(nComb,1);
res.RMSE = nan(nComb,1);
res.CORR = nan(nComb,1);

%% ===================== PARALELISMO + PROGRESSO =====================
useParallel = false;
try
    useParallel = license('test','Distrib_Computing_Toolbox') && exist('parpool','file')==2;
catch
    useParallel = false;
end

dq = [];
useDQ = false;
if useParallel
    try
        p = gcp('nocreate');
        if isempty(p), parpool; end
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(~)progress_update(nComb, PRINT_EVERY));
        useDQ = true;
        fprintf('Parallel ON (parfor)\n');
    catch
        useParallel = false;
    end
end
if ~useParallel
    fprintf('Parallel OFF (for)\n');
end

%% ===================== RODA GRID =====================
tAll = tic;

if useParallel
    parfor ic = 1:nComb
        % mapeia ic -> (idx0, idxMed)
        idx0   = ceil(ic / numel(MEDV));
        idxMed = ic - (idx0-1)*numel(MEDV);

        P = FIX;

        P.WAVELET        = GRID.WAVELET{comb0(idx0,1)};
        P.VOICES_PER_OCT = comb0(idx0,2);

        P.BP_ORDER = comb0(idx0,3);
        P.BP_WI    = comb0(idx0,4);
        P.BP_WF    = comb0(idx0,5);

        P.RIDGE_LAMBDA_HZ2 = comb0(idx0,6);
        P.HARM_RATIO       = comb0(idx0,7);
        P.MEDIAN_SMOOTH_SEC = MEDV(idxMed);

        out = run_one_pair(CSV_PATH, POLAR_HR_PATH, COL_T_MS, COL_HEART, gap_thr, P);

        res_wavelet_id(ic) = comb0(idx0,1);
        res_vpo(ic)        = P.VOICES_PER_OCT;
        res_bp_order(ic)   = P.BP_ORDER;
        res_bp_wi(ic)      = P.BP_WI;
        res_bp_wf(ic)      = P.BP_WF;
        res_lambda_hz2(ic) = P.RIDGE_LAMBDA_HZ2;
        res_harm_ratio(ic) = P.HARM_RATIO;
        res_med_sec(ic)    = P.MEDIAN_SMOOTH_SEC;

        res_MAE(ic)  = out.MAE;
        res_RMSE(ic) = out.RMSE;
        res_CORR(ic) = out.CORR;

        if useDQ, send(dq,1); end
    end

    % coleta do parfor (workaround: variáveis temporárias)
    res.wavelet_id = res_wavelet_id; %#ok<NASGU>
else
    count = 0;
    for idx0 = 1:n0
        for idxMed = 1:numel(MEDV)
            count = count + 1;

            P = FIX;
            P.WAVELET        = GRID.WAVELET{comb0(idx0,1)};
            P.VOICES_PER_OCT = comb0(idx0,2);

            P.BP_ORDER = comb0(idx0,3);
            P.BP_WI    = comb0(idx0,4);
            P.BP_WF    = comb0(idx0,5);

            P.RIDGE_LAMBDA_HZ2 = comb0(idx0,6);
            P.HARM_RATIO       = comb0(idx0,7);
            P.MEDIAN_SMOOTH_SEC = MEDV(idxMed);

            out = run_one_pair(CSV_PATH, POLAR_HR_PATH, COL_T_MS, COL_HEART, gap_thr, P);

            res.wavelet_id(count) = comb0(idx0,1);
            res.vpo(count)        = P.VOICES_PER_OCT;
            res.bp_order(count)   = P.BP_ORDER;
            res.bp_wi(count)      = P.BP_WI;
            res.bp_wf(count)      = P.BP_WF;
            res.lambda_hz2(count) = P.RIDGE_LAMBDA_HZ2;
            res.harm_ratio(count) = P.HARM_RATIO;
            res.med_sec(count)    = P.MEDIAN_SMOOTH_SEC;

            res.MAE(count)  = out.MAE;
            res.RMSE(count) = out.RMSE;
            res.CORR(count) = out.CORR;

            if mod(count, PRINT_EVERY)==0 || count==nComb
                progress_update(nComb, PRINT_EVERY, count, tAll);
            end
        end
    end
end

% Se rodou em parfor, precisamos puxar as temporárias do workspace do parfor:
if useParallel
    res.wavelet_id = res_wavelet_id;
    res.vpo        = res_vpo;
    res.bp_order   = res_bp_order;
    res.bp_wi      = res_bp_wi;
    res.bp_wf      = res_bp_wf;
    res.lambda_hz2 = res_lambda_hz2;
    res.harm_ratio = res_harm_ratio;
    res.med_sec    = res_med_sec;
    res.MAE        = res_MAE;
    res.RMSE       = res_RMSE;
    res.CORR       = res_CORR;
end

fprintf('\nGrid finalizado.\n');

%% ===================== BEST + TOP-10 =====================
valid = isfinite(res.RMSE);
if ~any(valid)
    error('Nenhuma combinação gerou métricas válidas (verifique CSV/Polar).');
end

idx = find(valid);
M = [res.RMSE(valid), res.MAE(valid), -res.CORR(valid)];
bestRow = lexicomin_row(M);
bestIdx = idx(bestRow);

wbest = GRID.WAVELET{res.wavelet_id(bestIdx)};

fprintf('\n===== MELHOR (menor RMSE; desempate MAE; depois maior CORR) =====\n');
fprintf('wavelet=%s | VPO=%d | BP=%d [%.2f %.2f] Hz | lambda_hz2=%.1f | harm=%.2f | med=%.1fs\n', ...
    wbest, res.vpo(bestIdx), res.bp_order(bestIdx), res.bp_wi(bestIdx), res.bp_wf(bestIdx), ...
    res.lambda_hz2(bestIdx), res.harm_ratio(bestIdx), res.med_sec(bestIdx));
fprintf('MAE=%.4f | RMSE=%.4f | CORR=%.4f\n', res.MAE(bestIdx), res.RMSE(bestIdx), res.CORR(bestIdx));

[~, ord] = sort(res.RMSE(valid), 'ascend');
ord = ord(1:min(10,numel(ord)));
topIdx = idx(ord);

fprintf('\n===== TOP-%d por RMSE =====\n', numel(topIdx));
for k=1:numel(topIdx)
    ii = topIdx(k);
    wk = GRID.WAVELET{res.wavelet_id(ii)};
    fprintf('#%02d RMSE=%.4f MAE=%.4f CORR=%.4f | w=%s vpo=%d | BP%d[%.2f %.2f] | lam=%.1f harm=%.2f med=%.1fs\n', ...
        k, res.RMSE(ii), res.MAE(ii), res.CORR(ii), wk, res.vpo(ii), ...
        res.bp_order(ii), res.bp_wi(ii), res.bp_wf(ii), res.lambda_hz2(ii), res.harm_ratio(ii), res.med_sec(ii));
end

%% ===================== PLOT DO MELHOR =====================
if FIX.PLOT_BEST
    P = FIX;
    P.WAVELET        = wbest;
    P.VOICES_PER_OCT = res.vpo(bestIdx);
    P.BP_ORDER       = res.bp_order(bestIdx);
    P.BP_WI          = res.bp_wi(bestIdx);
    P.BP_WF          = res.bp_wf(bestIdx);
    P.RIDGE_LAMBDA_HZ2= res.lambda_hz2(bestIdx);
    P.HARM_RATIO     = res.harm_ratio(bestIdx);
    P.MEDIAN_SMOOTH_SEC= res.med_sec(bestIdx);

    out = run_one_pair(CSV_PATH, POLAR_HR_PATH, COL_T_MS, COL_HEART, gap_thr, P, true); %#ok<NASGU>
end

fprintf('\nFim. (Sem save)\n');

%% =====================================================================
%% ===================== FUNÇÕES =====================
function out = run_one_pair(CSV_PATH, POLAR_HR_PATH, COL_T_MS, COL_HEART, gap_thr, P, doPlot)
    if nargin < 7, doPlot = false; end
    out = struct('MAE',NaN,'RMSE',NaN,'CORR',NaN);

    data = readmatrix(CSV_PATH);
    if size(data,2) < max(COL_T_MS, COL_HEART), return; end

    tpuro = data(:,COL_T_MS);
    x0    = data(:,COL_HEART);

    tpuro = tpuro(:); x0 = x0(:);
    okx = isfinite(x0);
    if ~all(okx)
        x0(~okx) = interp1(find(okx), x0(okx), find(~okx), 'linear', 'extrap');
    end

    t = (tpuro - tpuro(1))/1000;
    [t, ia] = unique(t, 'stable');
    x0 = x0(ia);

    nt = numel(t);
    if nt < 64, return; end

    dt_all = diff(t);
    break_ix = find(dt_all > gap_thr);

    if isempty(break_ix)
        seg_starts = 1; seg_ends = nt;
    else
        seg_starts = [1; break_ix+1];
        seg_ends   = [break_ix; nt];
    end

    hr_bpm = nan(nt,1);

    for s = 1:numel(seg_starts)
        idx = seg_starts(s):seg_ends(s);
        if numel(idx) < 64, continue; end

        t_seg = t(idx);
        x_seg = x0(idx);

        dt_seg = diff(t_seg);
        if any(dt_seg <= 0), continue; end
        fs = 1/mean(dt_seg);

        % DC
        x_seg = x_seg - mean(x_seg);

        % ====== FILTRO (NOVO) ======
        if ~isempty(P.BP_ORDER) && P.BP_ORDER > 0
            if ~(isfinite(P.BP_WI) && isfinite(P.BP_WF) && P.BP_WI > 0 && P.BP_WI < P.BP_WF && P.BP_WF < fs/2)
                continue;
            end
            [b,a] = butter(P.BP_ORDER, [P.BP_WI P.BP_WF]/(fs/2), 'bandpass');
            if any(~isfinite(x_seg)), continue; end
            x_seg = filtfilt(b,a,x_seg);
        end

        if any(~isfinite(x_seg)), continue; end

        % CWT compatível
        [wt,f] = get_cwt_compat(x_seg, fs, P.WAVELET, P.VOICES_PER_OCT);
        f = f(:);
        [f, ord] = sort(f,'ascend');
        wt = wt(ord,:);

        Pow = abs(wt).^2;
        if isempty(Pow) || size(Pow,2) ~= numel(idx), continue; end

        % evita DC (bin 1)
        Pow(1,:) = 0;

        Pnorm = Pow ./ (max(Pow,[],1) + eps);

        ridge_f = ridge_track_continuity_hz(Pnorm, f, P.RIDGE_LAMBDA_HZ2);
        ridge_f = harmonic_fix_half(Pow, f, ridge_f, P.HARM_RATIO);

        bpm = 60 * ridge_f(:);

        % suavização por segmento (mediana + movmean FIXO)
        wMed  = max(3, round(P.MEDIAN_SMOOTH_SEC * fs));
        wMean = max(3, round(P.MEAN_SMOOTH_SEC   * fs));

        bpm = movmedian(bpm, wMed, 'omitnan');
        if P.FILL_NAN_WITHIN_SEG
            bpm = fill_nan_by_interp(t_seg, bpm);
        end
        bpm = movmean(bpm, wMean, 'omitnan');

        hr_bpm(idx) = bpm;
    end

    if ~P.DO_POLAR_COMPARE || ~isfile(POLAR_HR_PATH)
        if doPlot
            figure('Color','w'); plot(t, hr_bpm, 'k', 'LineWidth', 1.4); grid on;
            xlabel('t (s)'); ylabel('HR (bpm)');
            title(sprintf('Radar HR | wavelet=%s vpo=%d | BP%d[%.2f %.2f] | lam=%.1f harm=%.2f med=%.1fs', ...
                P.WAVELET, P.VOICES_PER_OCT, P.BP_ORDER, P.BP_WI, P.BP_WF, P.RIDGE_LAMBDA_HZ2, P.HARM_RATIO, P.MEDIAN_SMOOTH_SEC));
        end
        return;
    end

    [t_polar, HR_polar] = read_txt_polar(POLAR_HR_PATH);
    [t_polar, ia] = unique(t_polar, 'stable');
    HR_polar = HR_polar(ia);

    okR = isfinite(t) & isfinite(hr_bpm);
    okP = isfinite(t_polar) & isfinite(HR_polar);
    if nnz(okR) < 2 || nnz(okP) < 2, return; end

    t_start = max(t(find(okR,1,'first')), t_polar(find(okP,1,'first')));
    t_end   = min(t(find(okR,1,'last')),  t_polar(find(okP,1,'last')));
    if ~(isfinite(t_start) && isfinite(t_end)) || (t_end <= t_start), return; end

    t_common = (t_start:P.DT_GRID:t_end)';
    HR_radar_grid = interp1(t(okR), hr_bpm(okR), t_common, 'linear');
    HR_polar_grid = interp1(t_polar(okP), HR_polar(okP), t_common, 'linear');

    mask_time = (t_common >= (t_start + P.IGNORE_START_SEC));
    HR_ref  = HR_polar_grid(mask_time);
    HR_test = HR_radar_grid(mask_time);

    ok = isfinite(HR_ref) & isfinite(HR_test);
    if nnz(ok) < 10, return; end

    err = HR_test(ok) - HR_ref(ok);

    out.MAE  = mean(abs(err), 'omitnan');
    out.RMSE = sqrt(mean(err.^2, 'omitnan'));
    out.CORR = corr(HR_test(ok), HR_ref(ok), 'Rows','complete');

    if doPlot
        figure('Color','w');
        plot(t_polar, HR_polar, 'm.-', 'LineWidth', 1, 'MarkerSize', 8); hold on;
        plot(t, hr_bpm, 'k', 'LineWidth', 1.4);
        grid on; xlabel('t (s)'); ylabel('HR (bpm)');
        legend('Polar','Radar','Location','best');
        title(sprintf('BEST | MAE=%.3f RMSE=%.3f CORR=%.3f | w=%s vpo=%d | BP%d[%.2f %.2f] | lam=%.1f harm=%.2f med=%.1fs', ...
            out.MAE, out.RMSE, out.CORR, P.WAVELET, P.VOICES_PER_OCT, P.BP_ORDER, P.BP_WI, P.BP_WF, P.RIDGE_LAMBDA_HZ2, P.HARM_RATIO, P.MEDIAN_SMOOTH_SEC));
    end
end

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

        try, [cfs,f] = cwt(fb, x); return; catch, end
        try, [cfs,f] = cwt(x, 'FilterBank', fb); return; catch, end
        try, [cfs,f] = fb.cwt(x); return; catch, end

        error('cwtfilterbank existe, mas não encontrei forma compatível de aplicar a CWT.');
    end

    try, [cfs,f] = cwt(x, fs); return; catch, end
    try, [cfs,f] = cwt(x, 1/fs); return; catch, end
    error('Sem suporte a CWT nesta instalação.');
end

function ridge_f = ridge_track_continuity_hz(Pnorm, f, lambda_hz2)
    [~,L] = size(Pnorm);
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
        p1 = P(i1,k); p2 = P(i2,k);
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
    if fid == -1, error('Não foi possível abrir o Polar: %s', p); end
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

function progress_update(nTotal, printEvery, count, tAll)
    if nargin < 3
        % usado via DataQueue
        persistent c t0
        if isempty(c), c = 0; t0 = tic; end
        c = c + 1;
        if mod(c, printEvery)==0 || c==nTotal
            dt = toc(t0);
            rate = c / max(dt, eps);
            eta  = (nTotal - c) / max(rate, eps);
            fprintf('Progresso: %d/%d (%.1f%%) | %.2f comb/s | ETA ~ %.1fs\n', ...
                c, nTotal, 100*c/nTotal, rate, eta);
        end
    else
        dt = toc(tAll);
        rate = count / max(dt, eps);
        eta  = (nTotal - count) / max(rate, eps);
        fprintf('Progresso: %d/%d (%.1f%%) | %.2f comb/s | ETA ~ %.1fs\n', ...
            count, nTotal, 100*count/nTotal, rate, eta);
    end
end

function row = lexicomin_row(M)
    [~, order] = sortrows(M, 1:size(M,2));
    row = order(1);
end
