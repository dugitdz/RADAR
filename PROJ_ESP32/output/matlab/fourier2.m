% ========================================================================
% PYTHON TRACKER (ALGORITHM) -> MATLAB (OFFLINE) + COMPARAÇÃO POLAR + PLOT
% - Sem CSV
% - Sem GUI realtime
% - Faz: gap reset + cheby2 bandpass causal + tracker FFT (prior+EMA)
% - Lê RADAR (csv) + POLAR (txt), compara e plota em ABAS (tabs)
% ========================================================================
clc; clear; close all;

%% ===================== PAIRS (RADAR + POLAR) =====================
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

% (opcional)
% pairs(4).name = "TESTE";
% pairs(4).radar_csv = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\Radar_teste.csv";
% pairs(4).polar_txt = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\Polar_teste.txt";

COL_T_MS  = 1;
COL_PHASE = 2;

%% ===================== CONFIG (espelho do Python) =====================
CFG = struct();

CFG.HR_COL = 1;   % 1=total 2=resp 3=heart
CFG.BR_COL = 2;

CFG.N      = 256;  % HR
CFG.hop    = 32;
CFG.N_BR   = 512;  % BR
CFG.hop_BR = 64;

CFG.FS_TARGET = 50.0;   % dataset -> resample p/ fs fixo
CFG.GAP_THR   = 0.6;
IMETH         = "linear";

% --- HR bandpass causal (cheby2, ordem FINAL par) ---
CFG.HR_FILT_ORD = 4;
CFG.HR_FILT_WI  = 0.595;
CFG.HR_FILT_WF  = 3.05;
CFG.HR_FILT_RS  = 30.50;

% --- BR bandpass causal (cheby2, ordem FINAL par) ---
CFG.BR_FILT_ORD = 4;
CFG.BR_FILT_WI  = 0.10;
CFG.BR_FILT_WF  = 0.50;
CFG.BR_FILT_RS  = 30.0;

% prior HR (Hz)
CFG.HR_PRIOR_INIT_BPM = 80.0;
CFG.HR_PRIOR_MIN_BPM  = 35.0;
CFG.HR_PRIOR_MAX_BPM  = 240.0;
CFG.HR_PRIOR_SIGMA_HZ = 0.40;
CFG.HR_PRIOR_ALPHA    = 0.05;
CFG.HR_EMA_TAU_SEC    = 4.0;

% prior BR (Hz)
CFG.BR_PRIOR_INIT_BRPM = 15.0;
CFG.BR_PRIOR_MIN_BRPM  = 4.0;
CFG.BR_PRIOR_MAX_BRPM  = 40.0;
CFG.BR_PRIOR_SIGMA_HZ  = 0.10;
CFG.BR_PRIOR_ALPHA     = 0.10;
CFG.BR_EMA_TAU_SEC     = 6.0;

CFG.F_MIN_HZ = 0.05;

MAXP = 2000;

% grade comum para métricas/plot (igual teu exemplo)
DT_GRID = 0.1;

%% ===================== filtro (NO SEU ESTILO: design_bp -> b,a) =====================
[b_hr, a_hr] = design_bp(CFG.FS_TARGET, "cheby2", CFG.HR_FILT_ORD, CFG.HR_FILT_WI, CFG.HR_FILT_WF, CFG.HR_FILT_RS);
[b_br, a_br] = design_bp(CFG.FS_TARGET, "cheby2", CFG.BR_FILT_ORD, CFG.BR_FILT_WI, CFG.BR_FILT_WF, CFG.BR_FILT_RS);

% estado inicial para filter(b,a,.,zi) (equivalente prático ao zi0*x_init)
zi_hr_unit = zi_for_constant_input(b_hr, a_hr, 1.0);
zi_br_unit = zi_for_constant_input(b_br, a_br, 1.0);

%% ===================== FIGURA COM ABAS =====================
fig = figure('Name','PY-Tracker vs Polar (tabs)','Color','w','Position',[60 60 1500 850]);
tg  = uitabgroup(fig);

tabs = gobjects(numel(pairs),1);
axs  = gobjects(numel(pairs),1);
for d=1:numel(pairs)
    tabs(d) = uitab(tg,'Title',char(pairs(d).name));
    axs(d)  = axes('Parent',tabs(d)); hold(axs(d),'on'); grid(axs(d),'on');
end

%% ===================== LOOP =====================
nDS = numel(pairs);
MAEds  = nan(nDS,1);
RMSEds = nan(nDS,1);
CORRds = nan(nDS,1);

fprintf("\n=== PY TRACKER (offline) + Polar ===\n");
fprintf("FS=%.2f | GAP=%.2fs | HR N=%d hop=%d | BR N=%d hop=%d\n", ...
    CFG.FS_TARGET, CFG.GAP_THR, CFG.N, CFG.hop, CFG.N_BR, CFG.hop_BR);
fprintf("HR BP cheby2 ord=%d wi=%.3f wf=%.3f | BR BP cheby2 ord=%d wi=%.3f wf=%.3f\n\n", ...
    CFG.HR_FILT_ORD, CFG.HR_FILT_WI, CFG.HR_FILT_WF, ...
    CFG.BR_FILT_ORD, CFG.BR_FILT_WI, CFG.BR_FILT_WF);

for d = 1:nDS
    name = string(pairs(d).name);

    % ---------- RADAR ----------
    A = try_readmatrix(pairs(d).radar_csv);
    if size(A,2) < 2
        error("[%s] CSV com poucas colunas: %s", name, pairs(d).radar_csv);
    end

    if size(A,2) >= 4
        t_ms  = A(:,1);
        total = A(:,2);
        resp  = A(:,3);
        heart = A(:,4);
    else
        t_ms  = A(:,COL_T_MS);
        ph    = A(:,COL_PHASE);
        total = ph; resp = ph; heart = ph;
    end

    t_ms  = double(t_ms(:));
    total = force_finite_vector(double(total(:)));
    resp  = force_finite_vector(double(resp(:)));
    heart = force_finite_vector(double(heart(:)));

    ok = isfinite(t_ms);
    t_ms  = t_ms(ok);
    total = total(ok);
    resp  = resp(ok);
    heart = heart(ok);

    if isempty(t_ms)
        error("[%s] radar vazio", name);
    end

    t0 = (t_ms - t_ms(1))/1000.0;
    [t0, iu] = unique(t0,'stable');
    total = total(iu); resp = resp(iu); heart = heart(iu);

    segs = segment_by_gaps(t0, CFG.GAP_THR);
    dt   = 1/CFG.FS_TARGET;

    % trackers
    trk_hr = make_tracker(CFG.FS_TARGET, CFG.N, CFG.hop, CFG.F_MIN_HZ, ...
        CFG.HR_PRIOR_INIT_BPM, CFG.HR_PRIOR_MIN_BPM, CFG.HR_PRIOR_MAX_BPM, ...
        CFG.HR_PRIOR_SIGMA_HZ, CFG.HR_PRIOR_ALPHA, CFG.HR_EMA_TAU_SEC, MAXP);

    trk_br = make_tracker(CFG.FS_TARGET, CFG.N_BR, CFG.hop_BR, CFG.F_MIN_HZ, ...
        CFG.BR_PRIOR_INIT_BRPM, CFG.BR_PRIOR_MIN_BRPM, CFG.BR_PRIOR_MAX_BRPM, ...
        CFG.BR_PRIOR_SIGMA_HZ, CFG.BR_PRIOR_ALPHA, CFG.BR_EMA_TAU_SEC, MAXP);

    t_hr = []; y_hr = [];   % saídas HR (bpm)
    % (BR fica calculada mas não é usada nas métricas aqui)
    % t_br = []; y_br = [];

    for s = 1:size(segs,1)
        i0s = segs(s,1); i1s = segs(s,2);

        ts  = t0(i0s:i1s);
        tot = total(i0s:i1s);
        rr  = resp(i0s:i1s);
        hh  = heart(i0s:i1s);

        [ts, iu2] = unique(ts,'stable');
        tot = tot(iu2); rr = rr(iu2); hh = hh(iu2);
        if numel(ts) < 8, continue; end

        % resample uniforme
        tnew = (ts(1):dt:ts(end))';
        tot_u = interp1(ts, tot, tnew, IMETH);
        rr_u  = interp1(ts, rr,  tnew, IMETH);
        hh_u  = interp1(ts, hh,  tnew, IMETH);

        ok2 = isfinite(tnew) & isfinite(tot_u) & isfinite(rr_u) & isfinite(hh_u);
        tnew  = tnew(ok2);
        tot_u = force_finite_vector(tot_u(ok2));
        rr_u  = force_finite_vector(rr_u(ok2));
        hh_u  = force_finite_vector(hh_u(ok2));

        if numel(tnew) < max(CFG.N, CFG.N_BR) + 2, continue; end

        % reset_segment (igual Python)
        x0_hr = pick_col(tot_u(1), rr_u(1), hh_u(1), CFG.HR_COL);
        x0_br = pick_col(tot_u(1), rr_u(1), hh_u(1), CFG.BR_COL);

        zi_hr = zi_hr_unit * x0_hr;
        zi_br = zi_br_unit * x0_br;

        trk_hr = tracker_reset(trk_hr);
        trk_br = tracker_reset(trk_br);

        for i = 1:numel(tnew)
            t_sec = tnew(i);

            x_hr_raw = pick_col(tot_u(i), rr_u(i), hh_u(i), CFG.HR_COL);
            x_br_raw = pick_col(tot_u(i), rr_u(i), hh_u(i), CFG.BR_COL);

            % filtro causal (realtime-like)
            [yhr, zi_hr] = filter(b_hr, a_hr, x_hr_raw, zi_hr);
            [ybr, zi_br] = filter(b_br, a_br, x_br_raw, zi_br);

            [trk_hr, hr_est] = tracker_push(trk_hr, t_sec, yhr);
            [trk_br, ~]      = tracker_push(trk_br, t_sec, ybr);

            if ~isnan(hr_est)
                t_hr(end+1,1) = trk_hr.t_out(trk_hr.out_len);
                y_hr(end+1,1) = trk_hr.y_out(trk_hr.out_len);
            end
        end
    end

    if isempty(t_hr)
        error("[%s] tracker HR saiu vazio", name);
    end

    % ---------- POLAR ----------
    [t_ref, HR_ref] = read_txt_polar_flex(pairs(d).polar_txt);
    if isempty(t_ref) || isempty(HR_ref)
        error("[%s] polar vazio/erro: %s", name, pairs(d).polar_txt);
    end

    % ---------- COMPARAÇÃO (grade comum) ----------
    t_start = max(min(t_hr), min(t_ref));
    t_end   = min(max(t_hr), max(t_ref));
    if ~(isfinite(t_start) && isfinite(t_end) && t_end > t_start)
        error("[%s] sem interseção temporal radar/polar", name);
    end

    t_common = (t_start:DT_GRID:t_end)';

    HR_radar_grid = interp1_nan(t_hr,  y_hr,   t_common, IMETH);
    HR_polar_grid = interp1_nan(t_ref, HR_ref, t_common, IMETH);

    [MAE, RMSE, CORR, Nuse] = metrics(HR_radar_grid, HR_polar_grid);

    MAEds(d)  = MAE;
    RMSEds(d) = RMSE;
    CORRds(d) = CORR;

    fprintf("[%s] N=%d | MAE=%.3f | RMSE=%.3f | CORR=%.4f\n", name, Nuse, MAE, RMSE, CORR);

    % ---------- PLOT ----------
    ax = axs(d);
    cla(ax);
    plot(ax, t_common, HR_polar_grid, 'LineWidth', 1.6);
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

function x = pick_col(total, resp, heart, col_idx)
    if col_idx == 1
        x = double(total);
    elseif col_idx == 2
        x = double(resp);
    else
        x = double(heart);
    end
end

% filtro no teu estilo: design_bp(fs, "cheby2", ord_final, wi, wf, rs)
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

% zi para filter(b,a,.,zi) ~ "zi0 * x_init"
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

% ---------- tracker (FFT + prior + EMA) ----------
function trk = make_tracker(fs, N, hop, f_min_hz, ...
    prior_init_per_min, prior_min_per_min, prior_max_per_min, ...
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

    trk.prior_min = double(prior_min_per_min)/60.0;
    trk.prior_max = double(prior_max_per_min)/60.0;
    trk.f_prev = min(max(double(prior_init_per_min)/60.0, trk.prior_min), trk.prior_max);

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

function trk = tracker_reset(trk)
    trk.len = int32(0);
    trk.cnt = int32(0);
    trk.ema_state = NaN;
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
    k0bin = (double(i0)-1) + dp;   % 0-based fracionário
    f_est = k0bin * trk.fs / double(trk.N);

    if ~isfinite(f_est) || f_est <= 0, return; end

    y_inst = 60.0 * f_est;

    if isnan(trk.ema_state)
        trk.ema_state = y_inst;
    else
        trk.ema_state = (1.0 - trk.ema_alpha)*trk.ema_state + trk.ema_alpha*y_inst;
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