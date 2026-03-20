clc; clear; close all;

pairs = struct([]);

pairs(1).name = "H2";
pairs(1).radar_csv = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases.csv";
pairs(1).polar_txt = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH2.txt";

pairs(2).name = "H3";
pairs(2).radar_csv = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv";
pairs(2).polar_txt = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH3.txt";

pairs(3).name = "H10";
pairs(3).radar_csv = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phase.csv";
pairs(3).polar_txt = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH10.txt";

pairs(4).name = "TESTE";
pairs(4).radar_csv = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\radar_proto.csv";
pairs(4).polar_txt = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\polar_proto.txt";

COL_T_MS  = 1;
COL_PHASE = 2;

FS_TARGET = 50.0;
GAP_THR   = 0.6;
DT_GRID   = 0.1;
F_MIN_HZ  = 0.05;
IMETH     = "linear";

CFG = struct();
CFG.N    = 512;
CFG.hop  = 64;

CFG.FILT_ORD = 4;
CFG.FILT_WI  = 0.595;
CFG.FILT_WF  = 3.05;
CFG.FILT_RS  = 30.50;

CFG.PRIOR_INIT_BPM = 80.0;
CFG.PRIOR_MIN_BPM  = 35.0;
CFG.PRIOR_MAX_BPM  = 240.0;
CFG.PRIOR_SIGMA_HZ = 0.40;
CFG.PRIOR_ALPHA    = 0.05;
CFG.RT_EMA_TAU_SEC = 4.0;

MAXP = 2000;

[b, a] = design_bp(FS_TARGET, "cheby2", CFG.FILT_ORD, CFG.FILT_WI, CFG.FILT_WF, CFG.FILT_RS);
zi_unit = zi_for_constant_input(b, a, 1.0);

fig = figure('Name','PY-Tracker vs Polar (NO reset gaps | mean-before-filter)','Color','w','Position',[60 60 1500 850]);
tg  = uitabgroup(fig);

tabs = gobjects(numel(pairs),1);
axs  = gobjects(numel(pairs),1);
for d=1:numel(pairs)
    tabs(d) = uitab(tg,'Title',char(pairs(d).name));
    axs(d)  = axes('Parent',tabs(d)); hold(axs(d),'on'); grid(axs(d),'on');
end

nDS = numel(pairs);
MAEds  = nan(nDS,1);
RMSEds = nan(nDS,1);
CORRds = nan(nDS,1);

fprintf("\n=== MATLAB (no reset gaps + mean-before-filter) ===\n");
fprintf("FS=%.2f | GAP=%.2fs | N=%d hop=%d | BP cheby2 ord=%d wi=%.3f wf=%.3f\n\n", ...
    FS_TARGET, GAP_THR, CFG.N, CFG.hop, CFG.FILT_ORD, CFG.FILT_WI, CFG.FILT_WF);

for d = 1:nDS
    name = string(pairs(d).name);

    A = try_readmatrix(pairs(d).radar_csv);
    if size(A,2) < 2
        error("[%s] CSV com poucas colunas: %s", name, pairs(d).radar_csv);
    end

    t_ms = double(A(:,1));
    xraw = double(A(:,2));

    t_ms = t_ms(:); xraw = xraw(:);
    m = isfinite(t_ms) & isfinite(xraw);
    t_ms = t_ms(m); xraw = xraw(m);

    if numel(t_ms) < CFG.N
        error("[%s] poucos pontos no radar (n=%d)", name, numel(t_ms));
    end

    t_sec = (t_ms - t_ms(1))/1000.0;
    [t_sec, iu] = unique(t_sec,'stable');
    xraw = xraw(iu);

    x0 = xraw - mean(xraw, 'omitnan');
    x0 = force_finite_vector(x0);

    segs = segment_by_gaps(t_sec, GAP_THR);

    [t_u, x_u] = resample_segments_irregular_fast(t_sec, x0, segs, FS_TARGET, IMETH);

    if numel(t_u) < CFG.N
        error("[%s] poucos pontos após resample (n=%d)", name, numel(t_u));
    end

    [t_ref, HR_ref] = read_txt_polar_flex(pairs(d).polar_txt);
    if isempty(t_ref) || isempty(HR_ref) || numel(t_ref) < 5
        error("[%s] Polar vazio/curto: %s", name, pairs(d).polar_txt);
    end

    trk = make_tracker(FS_TARGET, CFG.N, CFG.hop, F_MIN_HZ, ...
        CFG.PRIOR_INIT_BPM, CFG.PRIOR_MIN_BPM, CFG.PRIOR_MAX_BPM, ...
        CFG.PRIOR_SIGMA_HZ, CFG.PRIOR_ALPHA, CFG.RT_EMA_TAU_SEC, MAXP);

    zi = zi_unit * x_u(1);

    t_out = [];
    hr_out = [];

    for i = 1:numel(x_u)
        [y, zi] = filter(b, a, x_u(i), zi);
        [trk, hr_est] = tracker_push(trk, t_u(i), y);

        if ~isnan(hr_est)
            t_out(end+1,1)  = trk.t_out(trk.out_len);
            hr_out(end+1,1) = trk.y_out(trk.out_len);
        end
    end

    if numel(t_out) < 3
        error("[%s] tracker retornou muito pouco (n=%d)", name, numel(t_out));
    end

    t_start = max(min(t_out), min(t_ref));
    t_end   = min(max(t_out), max(t_ref));
    if ~(isfinite(t_start) && isfinite(t_end) && t_end > t_start)
        error("[%s] sem interseção temporal radar/polar", name);
    end

    t_common = (t_start:DT_GRID:t_end)';

    HR_radar_grid = interp1_nan(t_out, hr_out, t_common, IMETH);
    HR_polar_grid = interp1_nan(t_ref, HR_ref, t_common, IMETH);

    [MAE, RMSE, CORR, Nuse] = metrics(HR_radar_grid, HR_polar_grid);

    MAEds(d)  = MAE;
    RMSEds(d) = RMSE;
    CORRds(d) = CORR;

    fprintf("[%s] N=%d | MAE=%.3f | RMSE=%.3f | CORR=%.4f\n", name, Nuse, MAE, RMSE, CORR);

    ax = axs(d);
    cla(ax);
    plot(ax, t_common, HR_polar_grid, 'LineWidth', 1.6); hold(ax,'on');
    plot(ax, t_common, HR_radar_grid, 'LineWidth', 1.4);
    grid(ax,'on');
    xlabel(ax,'t (s)');
    ylabel(ax,'HR (bpm)');
    title(ax, sprintf('%s | MAE=%.2f RMSE=%.2f CORR=%.3f | N=%d', name, MAE, RMSE, CORR, Nuse));
    legend(ax,'POLAR','RADAR(PY-tracker)','Location','best');
end

fprintf("\n=== MÉDIA (%d datasets) ===\n", nDS);
fprintf("MAE=%.3f | RMSE=%.3f | CORR=%.4f\n", mean(MAEds,'omitnan'), mean(RMSEds,'omitnan'), mean(CORRds,'omitnan'));

tg.SelectedTab = tabs(1);

function A = try_readmatrix(p)
    try
        A = readmatrix(p);
    catch
        T = readtable(p);
        A = table2array(T);
    end
end

function x = force_finite_vector(x)
    x = double(x(:));
    ok = isfinite(x);
    if all(ok), return; end
    if nnz(ok) >= 2
        xi = find(ok);
        x(~ok) = interp1(xi, x(ok), find(~ok), 'linear', 'extrap');
    elseif nnz(ok) == 1
        x(~ok) = x(ok);
    else
        x(:) = 0; return;
    end
    x(~isfinite(x)) = 0;
end

function segments = segment_by_gaps(t, gap_thr)
    t = double(t(:));
    n = numel(t);
    if n < 2
        segments = [1 n];
        return;
    end
    brk = find(diff(t) > double(gap_thr));
    if isempty(brk)
        segments = [1 n];
    else
        segments = [[1; brk+1], [brk; n]];
    end
end

function [t_u, x_u] = resample_segments_irregular_fast(t, x, segments, fs_target, method)
    t = double(t(:)); x = double(x(:));
    dt = 1.0/double(fs_target);

    tCells = cell(size(segments,1),1);
    xCells = cell(size(segments,1),1);

    for k = 1:size(segments,1)
        i0 = segments(k,1); i1 = segments(k,2);
        ts = t(i0:i1); xs = x(i0:i1);

        [ts, iu] = unique(ts,'stable');
        xs = xs(iu);

        if numel(ts) < 4, continue; end

        tnew = (ts(1):dt:ts(end))';
        xnew = interp1(ts, xs, tnew, method);

        ok = isfinite(tnew) & isfinite(xnew);
        tCells{k} = tnew(ok);
        xCells{k} = xnew(ok);
    end

    if all(cellfun(@isempty,tCells))
        t_u = []; x_u = []; return;
    end

    t_u = vertcat(tCells{:});
    x_u = vertcat(xCells{:});

    [t_u, iu] = unique(t_u,'stable');
    x_u = x_u(iu);
end

function [b,a] = design_bp(fs, ftype, ord_final, wi, wf, rs)
    fs = double(fs); wi = double(wi); wf = double(wf);
    ord_final = double(ord_final);

    if mod(ord_final,2) ~= 0
        error("ord_final deve ser par");
    end
    n = ord_final/2;
    Wn = [wi wf]/(fs/2);

    switch lower(string(ftype))
        case "butter"
            [b,a] = butter(n, Wn, 'bandpass');
        case "cheby2"
            if nargin < 6 || isempty(rs), rs = 30; end
            [b,a] = cheby2(n, rs, Wn, 'bandpass');
        otherwise
            error("Filtro não suportado: %s", string(ftype));
    end
end

function zi = zi_for_constant_input(b,a, x0)
    b = b(:); a = a(:);
    na = numel(a); nb = numel(b);
    nf = max(na, nb);
    if nf < 2
        zi = [];
        return;
    end
    if nb < nf, b(end+1:nf) = 0; end
    if na < nf, a(end+1:nf) = 0; end

    nst = nf - 1;
    xpast = x0 * ones(nst,1);
    ypast = x0 * ones(nst,1);
    zi = filtic(b, a, ypast, xpast);
end

function trk = make_tracker(fs, N, hop, f_min_hz, ...
    prior_init_bpm, prior_min_bpm, prior_max_bpm, ...
    prior_sigma_hz, prior_alpha, ema_tau_sec, maxp)

    trk = struct();
    trk.fs  = double(fs);
    trk.N   = int32(N);
    trk.hop = int32(hop);
    trk.f_min_hz = double(f_min_hz);

    n = (0:double(N)-1)';
    trk.win = (0.21557895 ...
        - 0.41663158*cos(2*pi*n/double(N)) ...
        + 0.277263158*cos(4*pi*n/double(N)) ...
        - 0.083578947*cos(6*pi*n/double(N)) ...
        + 0.006947368*cos(8*pi*n/double(N)));

    trk.k_min = max(1, ceil(trk.f_min_hz * double(N) / trk.fs) - 1); % 0-based
    trk.f_bins = (0:floor(double(N)/2))' * (trk.fs/double(N));

    trk.prior_min = double(prior_min_bpm)/60.0;
    trk.prior_max = double(prior_max_bpm)/60.0;
    trk.f_prev = min(max(double(prior_init_bpm)/60.0, trk.prior_min), trk.prior_max);

    trk.prior_sigma_hz = double(prior_sigma_hz);
    trk.prior_alpha    = double(prior_alpha);

    if ~isempty(ema_tau_sec) && double(ema_tau_sec) > 0
        trk.ema_alpha = 1.0 - exp(-(double(hop)/trk.fs)/double(ema_tau_sec));
    else
        trk.ema_alpha = 1.0;
    end
    trk.ema_state = NaN;

    trk.buf  = zeros(double(N),1);
    trk.tbuf = zeros(double(N),1);
    trk.len  = int32(0);
    trk.cnt  = int32(0);

    trk.t_out   = zeros(double(maxp),1);
    trk.y_out   = zeros(double(maxp),1);
    trk.out_len = int32(0);
end

function [trk, y_per_min] = tracker_push(trk, t, x)
    y_per_min = NaN;

    if trk.len < trk.N
        trk.len = trk.len + 1;
        idx = double(trk.len);
        trk.buf(idx)  = double(x);
        trk.tbuf(idx) = double(t);
    else
        trk.buf(1:end-1)  = trk.buf(2:end);
        trk.tbuf(1:end-1) = trk.tbuf(2:end);
        trk.buf(end)      = double(x);
        trk.tbuf(end)     = double(t);
    end

    trk.cnt = trk.cnt + 1;
    if trk.len < trk.N, return; end
    if trk.cnt < trk.hop, return; end
    trk.cnt = int32(0);

    xw = trk.buf(:) .* trk.win(:);
    xw = xw - mean(xw);

    X  = fft(xw, double(trk.N));
    nh = floor(double(trk.N)/2) + 1;
    P  = abs(X(1:nh)).^2;

    k_srch_min0 = max(trk.k_min, ceil(trk.prior_min * double(trk.N) / trk.fs));
    k_srch_max0 = min(nh-1, floor(trk.prior_max * double(trk.N) / trk.fs));
    if k_srch_max0 <= k_srch_min0, return; end

    i_min = int32(k_srch_min0 + 1);
    i_max = int32(k_srch_max0 + 1);

    if isfinite(trk.prior_sigma_hz) && trk.prior_sigma_hz > 0
        W = exp(-0.5 * ((trk.f_bins - trk.f_prev) / trk.prior_sigma_hz).^2);
        tmp = P(:) .* W(:);
        [~, kk] = max(tmp(double(i_min):double(i_max)));
        i0 = double(i_min) - 1 + kk;
    else
        [~, kk] = max(P(double(trk.k_min+1):end));
        i0 = double(trk.k_min+1) - 1 + kk;
    end

    dp = parabolic_peak_logP(P, i0);
    k0bin = (double(i0)-1) + dp;
    f_est = k0bin * trk.fs / double(trk.N);

    if ~isfinite(f_est) || f_est <= 0, return; end

    hr_inst = 60.0 * f_est;

    if isnan(trk.ema_state)
        trk.ema_state = hr_inst;
    else
        trk.ema_state = (1.0 - trk.ema_alpha)*trk.ema_state + trk.ema_alpha*hr_inst;
    end

    f_clip = min(max(f_est, trk.prior_min), trk.prior_max);
    trk.f_prev = (1.0 - trk.prior_alpha)*trk.f_prev + trk.prior_alpha*f_clip;

    t_center = trk.tbuf(1) + 0.5*(double(trk.N)-1)/trk.fs;

    trk.out_len = trk.out_len + 1;
    if trk.out_len > int32(numel(trk.t_out))
        trk.t_out(1:end-1) = trk.t_out(2:end);
        trk.y_out(1:end-1) = trk.y_out(2:end);
        trk.out_len = int32(numel(trk.t_out));
        trk.t_out(end) = t_center;
        trk.y_out(end) = trk.ema_state;
    else
        trk.t_out(double(trk.out_len)) = t_center;
        trk.y_out(double(trk.out_len)) = trk.ema_state;
    end

    y_per_min = trk.ema_state;
end

function dp = parabolic_peak_logP(P, k)
    dp = 0.0;
    if k <= 1 || k >= numel(P), return; end
    v = log(P(k-1:k+1) + 1e-20);
    den = v(1) - 2*v(2) + v(3);
    if den == 0 || ~isfinite(den), return; end
    dp = 0.5*(v(1) - v(3))/den;
    dp = min(max(dp, -0.75), 0.75);
end

function [t_sec, HR] = read_txt_polar_flex(p)
    fid = fopen(p,'r');
    if fid == -1, error("Não abriu Polar: %s", p); end
    raw = textscan(fid,'%s','Delimiter','\n');
    raw = raw{1};
    fclose(fid);

    if numel(raw) < 2
        t_sec = []; HR = []; return;
    end

    lines = raw(2:end);
    ts_str = strings(numel(lines),1);
    HR = nan(numel(lines),1);

    for i = 1:numel(lines)
        L = string(lines{i});
        if contains(L,";"), parts = split(L,";"); else, parts = split(L,","); end
        if numel(parts) >= 2
            ts_str(i) = strtrim(parts(1));
            HR(i)     = str2double(strtrim(parts(2)));
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
    R = corr(yhat(m), yref(m), 'Rows','complete');
end