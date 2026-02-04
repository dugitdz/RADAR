% ========================================================================
% TOP1 CORR (3x POLAR: H2/H3/H10) — COM PARABÓLICA (logP)
% Plot em 3 ABAS (tabs) na MESMA janela, maior.
%
% Mantém o "teste sem 4.2" (sem máscara WI..WF no pico HR):
%   - pico HR escolhido no meio-espectro inteiro (pula DC)
% Volta a 4.4:
%   - parabólica em log(P) para refino sub-bin
% ========================================================================
clc; clear; close all;

%% ========================== PATHS (3 pares Polar) ==========================
pairs = struct([]);

pairs(1).name = "H2";
pairs(1).radar_csv = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phases.csv";
pairs(1).polar_txt = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\POLARH2.txt";

pairs(2).name = "H3";
pairs(2).radar_csv = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phases_raw.csv";
pairs(2).polar_txt = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\POLARH3.txt";

pairs(3).name = "H10";
pairs(3).radar_csv = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phase.csv";
pairs(3).polar_txt = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\POLARH10.txt";

%% ========================== FIXOS ==========================
FS_TARGET   = 25.0;
GAP_THR_SEC = 0.6;

DT_GRID   = 0.1;      % grade comum para comparar (sem suavizar)
IMETH     = "linear";

% score só pra print
ALPHA_CORR = 10;
BETA_MAE   = 0.2;

% low band (pré-BP)
LOW_F1_HZ = 0.10;
LOW_F2_HZ = 0.70;
LOW_NFFT  = 256;

% ridge continuidade (já estava no teu principal)
RIDGE_LAMBDA_HZ2 = 250;

%% ========================== TOP1 (da tua lista de CORR) ==========================
FILT_TYPE  = "cheby2";
ORD_FINAL  = 4;      % ordem FINAL (par)
WI         = 0.8;    % Hz (corte do bandpass)
WF         = 2.0;    % Hz

WIN_SPEC   = "rect";
N          = 32;
HOP        = 16;
NFFT       = 32;

SNR_THR_DB   = 6;
SNR_RANGE_DB = 9;
M_MAX        = 8;
TOL_HZ       = 0.08;
TOPK_HR      = 12;
BETA         = 3;

%% ========================== PREP: filtro / janela ==========================
[b_bp, a_bp] = design_bp(FS_TARGET, FILT_TYPE, ORD_FINAL, WI, WF);
win = make_window(N, WIN_SPEC);

% low_f2_eff (abaixo do WI)
low_f2_eff = min(LOW_F2_HZ, 0.90*WI);
if low_f2_eff <= LOW_F1_HZ
    error("LOW band inválida: low_f2_eff <= LOW_F1_HZ (WI=%.2f).", WI);
end

%% ========================== FIGURA COM 3 ABAS ==========================
fig = figure('Name','HR Radar vs Polar (3 abas)','Color','w', ...
             'Position',[60 60 1500 850]);
tg = uitabgroup(fig);

tabs = gobjects(numel(pairs),1);
axs  = gobjects(numel(pairs),1);
for d = 1:numel(pairs)
    tabs(d) = uitab(tg, 'Title', char(pairs(d).name));
    axs(d)  = axes('Parent', tabs(d));
end

%% ========================== RODA NOS 3 DATASETS ==========================
nDS = numel(pairs);
MAEds  = nan(nDS,1);
RMSEds = nan(nDS,1);
CORRds = nan(nDS,1);
SCOREds= nan(nDS,1);

fprintf("=== TOP1 CORR RUN (com parabólica; sem banda HR) ===\n");
fprintf("BP: %s ord=%d WI=%.2f WF=%.2f | WIN=%s | N=%d hop=%d NFFT=%d\n", ...
    FILT_TYPE, ORD_FINAL, WI, WF, WIN_SPEC, N, HOP, NFFT);
fprintf("LOW: [%.2f..%.2f]Hz (eff=%.2f) NFFT=%d | SNRthr=%g range=%g | M=%d tol=%.3f | TOPK=%d | beta=%g | ridge=%g\n\n", ...
    LOW_F1_HZ, LOW_F2_HZ, low_f2_eff, LOW_NFFT, SNR_THR_DB, SNR_RANGE_DB, M_MAX, TOL_HZ, TOPK_HR, BETA, RIDGE_LAMBDA_HZ2);

for d = 1:nDS
    name = string(pairs(d).name);
    fprintf("[%s] lendo radar...\n", name);

    % ----- RADAR -----
    A = try_readmatrix(pairs(d).radar_csv);
    if size(A,2) < 2
        error("CSV com poucas colunas: %s", pairs(d).radar_csv);
    end

    tpuro = A(:,1);
    xraw  = A(:,2);

    t0 = (tpuro - tpuro(1))/1000;  % s irregular
    t0 = t0(:); xraw = xraw(:);

    m = isfinite(t0) & isfinite(xraw);
    t0 = t0(m); xraw = xraw(m);

    x0 = prep_nowrap(xraw); % sem wrap

    seg0 = segment_by_gaps(t0, GAP_THR_SEC);
    [t_u, x_u] = resample_segments_irregular_fast(t0, x0, seg0, FS_TARGET, IMETH);

    if numel(t_u) < N
        error("[%s] poucos pontos após resample (n=%d < N=%d).", name, numel(t_u), N);
    end

    % pós-BP
    x_post = filtfilt_safe(b_bp, a_bp, x_u);

    % ----- POLAR -----
    fprintf("[%s] lendo polar...\n", name);
    [t_ref, HR_ref] = read_txt_polar_flex(pairs(d).polar_txt);
    if isempty(t_ref) || isempty(HR_ref)
        error("[%s] Polar vazio/erro: %s", name, pairs(d).polar_txt);
    end

    % ----- TRACK (sem banda HR, COM parabólica) -----
    fprintf("[%s] tracking...\n", name);
    [t_frames, f_hr] = rfft_track_dual_subharm_ridge_noband_parab( ...
        t_u, x_u, x_post, FS_TARGET, ...
        N, HOP, NFFT, win, ...
        LOW_F1_HZ, low_f2_eff, LOW_NFFT, ...
        SNR_THR_DB, SNR_RANGE_DB, ...
        M_MAX, TOL_HZ, TOPK_HR, BETA, ...
        RIDGE_LAMBDA_HZ2);

    if isempty(t_frames) || isempty(f_hr)
        error("[%s] tracker retornou vazio.", name);
    end

    HR_frames = 60*f_hr;

    % ----- INTERSEÇÃO / MÉTRICAS -----
    t_start = max(min(t_frames), min(t_ref));
    t_end   = min(max(t_frames), max(t_ref));
    if ~(isfinite(t_start) && isfinite(t_end) && t_end > t_start)
        error("[%s] não houve interseção temporal válida radar/polar.", name);
    end

    t_common = (t_start:DT_GRID:t_end)';

    HR_radar_grid = interp1_nan(t_frames, HR_frames, t_common, IMETH);
    HR_polar_grid = interp1_nan(t_ref,    HR_ref,   t_common, IMETH);

    [mae, rmse, rr, nuse] = metrics(HR_radar_grid, HR_polar_grid);
    if ~(isfinite(mae) && isfinite(rmse) && isfinite(rr) && nuse >= 3)
        error("[%s] métricas inválidas.", name);
    end

    MAEds(d)  = mae;
    RMSEds(d) = rmse;
    CORRds(d) = rr;
    SCOREds(d)= rmse + BETA_MAE*mae - ALPHA_CORR*rr;

    fprintf("[%s] n=%d | MAE=%.2f | RMSE=%.2f | CORR=%.3f | SCORE=%.3f\n\n", ...
        name, nuse, mae, rmse, rr, SCOREds(d));

    % ----- PLOT NA ABA -----
    ax = axs(d);
    cla(ax);
    plot(ax, t_common, HR_polar_grid, 'LineWidth', 1.5); hold(ax,'on');
    plot(ax, t_common, HR_radar_grid, 'LineWidth', 1.5);
    grid(ax,'on');
    xlabel(ax, "t (s)");
    ylabel(ax, "HR (bpm)");
    title(ax, sprintf("%s | MAE=%.2f RMSE=%.2f CORR=%.3f", name, mae, rmse, rr));
    legend(ax, "Polar","Radar","Location","best");
end

fprintf("=== MÉDIA (3 datasets) ===\n");
fprintf("MAE=%.3f | RMSE=%.3f | CORR=%.4f | SCORE=%.3f\n", ...
    mean(MAEds,'omitnan'), mean(RMSEds,'omitnan'), mean(CORRds,'omitnan'), mean(SCOREds,'omitnan'));

% abre na 1ª aba
tg.SelectedTab = tabs(1);

%% ========================================================================
%                               FUNÇÕES
% ========================================================================

function A = try_readmatrix(p)
    try
        A = readmatrix(p);
    catch
        T = readtable(p);
        A = table2array(T);
    end
end

function x = prep_nowrap(x)
    x = double(x(:));
    ok = isfinite(x);
    if ~all(ok)
        x(~ok) = interp1(find(ok), x(ok), find(~ok), 'linear', 'extrap');
    end
    x = x - mean(x,'omitnan');
end

function segments = segment_by_gaps(t, gap_thr)
    t = t(:);
    n = numel(t);
    if n < 2
        segments = [1 n];
        return;
    end
    brk = find(diff(t) > gap_thr);
    if isempty(brk)
        segments = [1 n];
    else
        segments = [[1; brk+1], [brk; n]];
    end
end

function [t_u, x_u] = resample_segments_irregular_fast(t, x, segments, fs_target, method)
    t = t(:); x = x(:);
    dt = 1/fs_target;

    tCells = cell(size(segments,1),1);
    xCells = cell(size(segments,1),1);

    for k = 1:size(segments,1)
        i0 = segments(k,1); i1 = segments(k,2);
        ts = t(i0:i1); xs = x(i0:i1);
        [ts, iu] = unique(ts,'stable'); xs = xs(iu);
        if numel(ts) < 4, continue; end

        tnew = (ts(1):dt:ts(end))';
        xnew = interp1(ts, xs, tnew, method);

        ok = isfinite(tnew) & isfinite(xnew);
        tCells{k} = tnew(ok);
        xCells{k} = xnew(ok);
    end

    t_u = vertcat(tCells{:});
    x_u = vertcat(xCells{:});

    [t_u, iu] = unique(t_u,'stable');
    x_u = x_u(iu);
end

function [t_sec, HR] = read_txt_polar_flex(p)
    fid = fopen(p,'r');
    if fid == -1, error("Não abriu Polar: %s", p); end

    raw = textscan(fid,'%s','Delimiter','\n');
    raw = raw{1};
    fclose(fid);

    if numel(raw) < 2
        t_sec = []; HR = [];
        return;
    end

    lines = raw(2:end);
    ts_str = strings(numel(lines),1);
    HR = nan(numel(lines),1);

    for i = 1:numel(lines)
        L = string(lines{i});
        if contains(L,";")
            parts = split(L,";");
        else
            parts = split(L,",");
        end
        if numel(parts) >= 2
            ts_str(i) = strtrim(parts(1));
            HR(i) = str2double(strtrim(parts(2)));
        end
    end

    try
        t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
    catch
        t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    end

    t_sec = seconds(t_dt - t_dt(1));
    t_sec = t_sec(:); HR = HR(:);

    m = isfinite(t_sec) & isfinite(HR);
    t_sec = t_sec(m); HR = HR(m);

    [t_sec, iu] = unique(t_sec,'stable');
    HR = HR(iu);
end

function [b,a] = design_bp(fs, ftype, ord_final, wi, wf)
    if mod(ord_final,2) ~= 0, error("ord_final deve ser par"); end
    n = ord_final/2;

    switch lower(string(ftype))
        case "butter"
            [b,a] = butter(n, [wi wf]/(fs/2), 'bandpass');
        case "cheby2"
            Rs = 30;
            [b,a] = cheby2(n, Rs, [wi wf]/(fs/2), 'bandpass');
        otherwise
            error("Filtro não suportado: %s", string(ftype));
    end
end

function y = filtfilt_safe(b,a,x)
    x = x(:);
    if numel(x) < 16
        y = x; return;
    end
    padlen = 3*(max(length(a), length(b)) - 1);
    if numel(x) <= padlen
        y = x; return;
    end
    y = filtfilt(b,a,x);
end

function w = make_window(N, spec)
    spec = lower(string(spec));
    parts = split(spec, ":");
    name = parts(1);

    hasParam = numel(parts) > 1;
    if hasParam
        p = str2double(parts(2));
        if ~isfinite(p), hasParam = false; end
    end

    switch name
        case "rect",     w = rectwin(N);
        case "hann",     w = hann(N,'periodic');
        case "hamming",  w = hamming(N,'periodic');
        case "blackman", w = blackman(N,'periodic');
        case "blackmanharris", w = blackmanharris(N,'periodic');
        case "flattop",  w = flattopwin(N,'periodic');
        case "nuttall",  w = nuttallwin(N,'periodic');
        case "kaiser"
            if ~hasParam, p = 8; end
            w = kaiser(N,p);
        case "tukey"
            if ~hasParam, p = 0.5; end
            w = tukeywin(N,p);
        case "cheb"
            if ~hasParam, p = 80; end
            w = chebwin(N,p);
        otherwise
            w = hann(N,'periodic');
    end
    w = w(:);
end

function [t_frames, f_hr] = rfft_track_dual_subharm_ridge_noband_parab( ...
    tu, x_pre, x_post, fs, ...
    nperseg, hop, nfft, win, ...
    low_f1, low_f2, low_nfft, ...
    snr_thr_db, snr_range_db, ...
    m_max, tol_hz, topk, beta, ...
    ridge_lambda_hz2)

    tu = tu(:); x_pre = x_pre(:); x_post = x_post(:);
    n = numel(x_post);

    if n < nperseg || numel(tu)~=n || numel(x_pre)~=n
        t_frames = []; f_hr = []; return;
    end

    nfft = max(nperseg, round(nfft));
    half = floor(nfft/2) + 1;

    % >>> sem banda HR (teste): usa tudo, pula DC
    k_min = 2;
    k_max = half;

    % low FFT (só se beta>0)
    hasLow = (beta > 0);
    if hasLow
        low_nfft = max(low_nfft, nperseg);
        halfL = floor(low_nfft/2) + 1;
        fbinL = (0:halfL-1)' * (fs/low_nfft);
        idxLow = find(fbinL >= low_f1 & fbinL <= low_f2);
        hasLow = ~isempty(idxLow);
    else
        idxLow = [];
    end

    t_frames = [];
    f_hr     = [];
    f_prev   = NaN;

    idx0 = 1;
    while (idx0 + nperseg - 1) <= n

        xp = x_pre(idx0:idx0+nperseg-1).*win;  xp = xp - mean(xp,'omitnan');
        xo = x_post(idx0:idx0+nperseg-1).*win; xo = xo - mean(xo,'omitnan');

        % -------- low (pré-BP) --------
        f_low = NaN; w_low = 0;
        if hasLow
            XL = fft(xp, low_nfft);
            PL = abs(XL(1:floor(low_nfft/2)+1)).^2;

            bandL = PL(idxLow);
            [pmaxL, relL] = max(bandL);
            kL = idxLow(1) + relL - 1;

            deltaL = parabolic_peak_interp_logP(PL, kL);
            f_low  = (kL - 1 + deltaL) * fs / low_nfft;

            pmedL  = median(bandL);
            snr_db = 10*log10((pmaxL+1e-30)/(pmedL+1e-30));

            w_low  = (snr_db - snr_thr_db)/snr_range_db;
            if ~isfinite(w_low), w_low = 0; end
            w_low  = max(0, min(1, w_low));
        end

        % -------- HR (pós-BP) --------
        X = fft(xo, nfft);
        P = abs(X(1:half)).^2;

        bandH = P(k_min:k_max);
        if isempty(bandH) || ~any(isfinite(bandH))
            idx0 = idx0 + hop;
            continue;
        end

        [~, ord] = sort(bandH, 'descend', 'MissingPlacement','last');
        ord = ord(1:min(topk, numel(ord)));
        cand_k = (k_min - 1) + ord(:);

        bestScore = -Inf;
        bestK = cand_k(1);

        for ii = 1:numel(cand_k)
            k0 = cand_k(ii);

            delta = parabolic_peak_interp_logP(P, k0);
            f_cand = (k0 - 1 + delta) * fs / nfft;

            base = log(P(k0) + 1e-30);

            pen = 0;
            if beta > 0 && w_low > 0 && isfinite(f_low)
                dmin = Inf;
                for m = 2:m_max
                    d = abs(f_cand - m*f_low);
                    if d < dmin, dmin = d; end
                end
                pen = exp(-(dmin/tol_hz)^2);
            end

            cont = 0;
            if isfinite(f_prev)
                cont = ridge_lambda_hz2 * (f_cand - f_prev).^2;
            end

            score = base - beta*w_low*pen - cont;

            if score > bestScore
                bestScore = score;
                bestK = k0;
            end
        end

        delta = parabolic_peak_interp_logP(P, bestK);
        f_est = (bestK - 1 + delta) * fs / nfft;

        if isfinite(f_est)
            f_prev = f_est;
        end

        t_center = tu(idx0) + 0.5*(nperseg-1)/fs;
        t_frames(end+1,1) = t_center; %#ok<AGROW>
        f_hr(end+1,1)     = f_est;    %#ok<AGROW>

        idx0 = idx0 + hop;
    end
end

function delta = parabolic_peak_interp_logP(P, k)
    if k <= 1 || k >= numel(P)
        delta = 0; return;
    end
    la = log(P(k-1) + 1e-30);
    lb = log(P(k)   + 1e-30);
    lc = log(P(k+1) + 1e-30);
    den = (la - 2*lb + lc);
    if den == 0 || ~isfinite(den)
        delta = 0; return;
    end
    delta = 0.5*(la - lc)/den;
    if ~isfinite(delta), delta = 0; end
    delta = max(-0.75, min(0.75, delta));
end

function y = interp1_nan(t, x, tq, method)
    t = t(:); x = x(:); tq = tq(:);
    ok = isfinite(t) & isfinite(x);
    y = nan(size(tq));
    if sum(ok) < 2, return; end
    try
        y = interp1(t(ok), x(ok), tq, method, nan);
    catch
        y = interp1(t(ok), x(ok), tq, 'linear', nan);
    end
end

function [MAE, RMSE, R, n] = metrics(yhat, yref)
    m = isfinite(yhat) & isfinite(yref);
    n = sum(m);
    if n < 3
        MAE = nan; RMSE = nan; R = nan; return;
    end
    e = yhat(m) - yref(m);
    MAE  = mean(abs(e));
    RMSE = sqrt(mean(e.^2));
    R = corr(yhat(m), yref(m));
end
