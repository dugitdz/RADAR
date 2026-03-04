clc; clear; close all;

% ===================== DATASETS (INPUT) =====================
BASE = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",    "radar", BASE+"phases.csv",      "polar", BASE+"POLARH2.txt")
    struct("name","H3",    "radar", BASE+"phases_raw.csv",  "polar", BASE+"POLARH3.txt")
    struct("name","H10",   "radar", BASE+"phase.csv",       "polar", BASE+"POLARH10.txt")
    struct("name","TESTE", "radar", BASE+"Radar_teste.csv", "polar", BASE+"Polar_teste.txt")
];

COL_T_MS  = 1;
COL_PHASE = 2;

% ===================== CSV OUTPUT =====================
CSV_DIR         = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\HR_BR";
CSV_ON          = true;
CSV_FLUSH_EVERY = 50;

% ===================== OFFLINE "REALTIME-LIKE" CONFIG =====================
CFG = struct();

% quais canais usar (1=total,2=resp,3=heart)
CFG.HR_COL = 1;
CFG.BR_COL = 2;

% N/hop separados
CFG.N      = 256;   % HR
CFG.hop    = 32;    % HR
CFG.N_BR   = 512;   % BR
CFG.hop_BR = 64;    % BR

% reamostragem via interp1 para FS_TARGET
CFG.FS_TARGET = 50.0;
CFG.GAP_THR   = 0.6;
IMETH         = 'linear';

% --- HR bandpass causal (cheby2, ordem PAR) ---
CFG.HR_FILT_ORD = 4;
CFG.HR_FILT_WI  = 0.595;
CFG.HR_FILT_WF  = 3.05;
CFG.HR_FILT_RS  = 30.50;

% --- BR bandpass causal (cheby2, ordem PAR) ---
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

% FFT min freq (evita bins muito baixos)
CFG.F_MIN_HZ = 0.05;

MAXP = 2000;

% ===================== BASIC VALIDATION =====================
if CFG.N < 32 || CFG.hop < 1 || CFG.hop > CFG.N
    error("CFG inválido (HR): verifique N/hop");
end
if CFG.N_BR < 32 || CFG.hop_BR < 1 || CFG.hop_BR > CFG.N_BR
    error("CFG inválido (BR): verifique N_BR/hop_BR");
end
if ~ismember(CFG.HR_COL,[1 2 3]) || ~ismember(CFG.BR_COL,[1 2 3])
    error("CFG inválido: HR_COL/BR_COL devem ser 1,2,3");
end

% ===================== FILTER DESIGN =====================
[sos_hr, zi_hr0, err] = design_sos_cheby2(CFG.FS_TARGET, CFG.HR_FILT_ORD, CFG.HR_FILT_WI, CFG.HR_FILT_WF, CFG.HR_FILT_RS);
if isempty(sos_hr), error("HR filter error: %s", err); end
[sos_br, zi_br0, err] = design_sos_cheby2(CFG.FS_TARGET, CFG.BR_FILT_ORD, CFG.BR_FILT_WI, CFG.BR_FILT_WF, CFG.BR_FILT_RS);
if isempty(sos_br), error("BR filter error: %s", err); end

% ===================== GUI (SIMPLE) =====================
fig = figure('Color','w','Name','Offline tracker (MATLAB) | dataset input');
tl  = tiledlayout(fig,3,1,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl,1); grid(ax1,'on'); hold(ax1,'on'); ylabel(ax1,'HR (bpm)');
ax2 = nexttile(tl,2); grid(ax2,'on'); hold(ax2,'on'); ylabel(ax2,'BR (brpm)');
ax3 = nexttile(tl,3); grid(ax3,'on'); hold(ax3,'on'); ylabel(ax3,'signals'); xlabel(ax3,'t (s)');
legend(ax3, {'sig_hr_raw','sig_br_raw','sig_hr_filt','sig_br_filt'}, 'Location','best');

fprintf("\n=== Offline tracker (MATLAB) ===\n");
fprintf("FS_TARGET=%.2f | GAP_THR=%.2fs\n", CFG.FS_TARGET, CFG.GAP_THR);
fprintf("HR: N=%d hop=%d | BP=[%.3f %.3f] Hz\n", CFG.N, CFG.hop, CFG.HR_FILT_WI, CFG.HR_FILT_WF);
fprintf("BR: N=%d hop=%d | BP=[%.3f %.3f] Hz\n", CFG.N_BR, CFG.hop_BR, CFG.BR_FILT_WI, CFG.BR_FILT_WF);
fprintf("=================================\n");

for d = 1:numel(DATASETS)
    ds = DATASETS(d);

    fprintf("\n---- %s ----\n", ds.name);

    A = readmatrix(ds.radar);

    % Aceita 2 colunas (t_ms, phase) OU 4 colunas (t_ms,total,resp,heart)
    if size(A,2) >= 4
        t_ms  = A(:,1);
        total = A(:,2);
        resp  = A(:,3);
        heart = A(:,4);
    else
        t_ms  = A(:,COL_T_MS);
        ph    = A(:,COL_PHASE);
        total = ph;
        resp  = ph;
        heart = ph;
        fprintf("Radar tem %d colunas -> usando phase como total/resp/heart (compat).\n", size(A,2));
    end

    t_ms  = double(t_ms(:));
    total = double(total(:));
    resp  = double(resp(:));
    heart = double(heart(:));

    ok = isfinite(t_ms) & isfinite(total) & isfinite(resp) & isfinite(heart);
    t_ms  = t_ms(ok);
    total = total(ok);
    resp  = resp(ok);
    heart = heart(ok);

    if isempty(t_ms)
        fprintf("%s | radar vazio\n", ds.name);
        continue;
    end

    % tempo em segundos a partir do primeiro timestamp
    t0 = (t_ms - t_ms(1))/1000.0;
    [t0, iu] = unique(t0, 'stable');
    t_ms  = t_ms(iu);
    total = total(iu);
    resp  = resp(iu);
    heart = heart(iu);

    % interp para FS_TARGET
    dt  = 1.0 / CFG.FS_TARGET;
    t_u = (t0(1):dt:t0(end))';

    total_u = interp1(t0, total, t_u, IMETH);
    resp_u  = interp1(t0, resp,  t_u, IMETH);
    heart_u = interp1(t0, heart, t_u, IMETH);

    total_u = force_finite_vector(total_u);
    resp_u  = force_finite_vector(resp_u);
    heart_u = force_finite_vector(heart_u);

    % segmentação por gaps (no timeline já uniformizado)
    segs = segment_by_gaps(t_u, CFG.GAP_THR);

    % trackers (reiniciados a cada segmento)
    trk_hr = make_tracker(CFG.FS_TARGET, CFG.N, CFG.hop, CFG.F_MIN_HZ, ...
        CFG.HR_PRIOR_INIT_BPM, CFG.HR_PRIOR_MIN_BPM, CFG.HR_PRIOR_MAX_BPM, ...
        CFG.HR_PRIOR_SIGMA_HZ, CFG.HR_PRIOR_ALPHA, CFG.HR_EMA_TAU_SEC, MAXP);

    trk_br = make_tracker(CFG.FS_TARGET, CFG.N_BR, CFG.hop_BR, CFG.F_MIN_HZ, ...
        CFG.BR_PRIOR_INIT_BRPM, CFG.BR_PRIOR_MIN_BRPM, CFG.BR_PRIOR_MAX_BRPM, ...
        CFG.BR_PRIOR_SIGMA_HZ, CFG.BR_PRIOR_ALPHA, CFG.BR_EMA_TAU_SEC, MAXP);

    last_hr = NaN;
    last_br = NaN;

    % CSV writer (arquivo por dataset)
    if CSV_ON
        if ~exist(CSV_DIR,'dir'), mkdir(CSV_DIR); end
        fname = "hr_br_" + ds.name + "_" + string(datetime("now","Format","yyyyMMdd_HHmmss")) + ".csv";
        csv_path = fullfile(CSV_DIR, fname);
        fid = fopen(csv_path, 'w');
        fprintf(fid, "timestamp,total,resp,heart,BR,HR\n");
        flush_cnt = 0;
        fprintf("CSV: %s\n", csv_path);
    else
        fid = -1;
    end

    % buffers para plot
    all_t_raw  = [];
    all_hr_raw = [];
    all_br_raw = [];
    all_hr_f   = [];
    all_br_f   = [];

    all_t_hr = [];
    all_hr   = [];
    all_t_br = [];
    all_br   = [];

    % filtros por segmento
    zi_hr = zi_hr0 .* 0;
    zi_br = zi_br0 .* 0;

    for s = 1:size(segs,1)
        i0 = segs(s,1); i1 = segs(s,2);

        % init de estado (tipo "reset_segment")
        x0_hr = pick_col(total_u(i0), resp_u(i0), heart_u(i0), CFG.HR_COL);
        x0_br = pick_col(total_u(i0), resp_u(i0), heart_u(i0), CFG.BR_COL);

        zi_hr  = zi_hr0 .* x0_hr;
        zi_br  = zi_br0 .* x0_br;
        trk_hr = tracker_reset(trk_hr);
        trk_br = tracker_reset(trk_br);

        for i = i0:i1
            t_sec = t_u(i);

            x_hr_raw = pick_col(total_u(i), resp_u(i), heart_u(i), CFG.HR_COL);
            x_br_raw = pick_col(total_u(i), resp_u(i), heart_u(i), CFG.BR_COL);

            [y_hr, zi_hr] = sosfilt(sos_hr, x_hr_raw, zi_hr);
            [y_br, zi_br] = sosfilt(sos_br, x_br_raw, zi_br);
            y_hr = double(y_hr);
            y_br = double(y_br);

            [trk_hr, hr_est] = tracker_push(trk_hr, t_sec, y_hr);
            [trk_br, br_est] = tracker_push(trk_br, t_sec, y_br);

            all_t_raw(end+1,1)  = t_sec;
            all_hr_raw(end+1,1) = x_hr_raw;
            all_br_raw(end+1,1) = x_br_raw;
            all_hr_f(end+1,1)   = y_hr;
            all_br_f(end+1,1)   = y_br;

            if ~isnan(hr_est)
                last_hr = hr_est;
                all_t_hr(end+1,1) = trk_hr.t_out(trk_hr.out_len);
                all_hr(end+1,1)   = trk_hr.y_out(trk_hr.out_len);
            end
            if ~isnan(br_est)
                last_br = br_est;
                all_t_br(end+1,1) = trk_br.t_out(trk_br.out_len);
                all_br(end+1,1)   = trk_br.y_out(trk_br.out_len);
            end

            if CSV_ON
                fprintf(fid, "%.3f,%.9g,%.9g,%.9g,%.9g,%.9g\n", ...
                    t_ms_from_tsec(t_sec, t_ms(1)), ...
                    pick_col(total_u(i), resp_u(i), heart_u(i), 1), ...
                    pick_col(total_u(i), resp_u(i), heart_u(i), 2), ...
                    pick_col(total_u(i), resp_u(i), heart_u(i), 3), ...
                    last_br, last_hr);

                flush_cnt = flush_cnt + 1;
                if CSV_FLUSH_EVERY > 0 && mod(flush_cnt, CSV_FLUSH_EVERY) == 0
                    fflush(fid);
                end
            end
        end
    end

    if CSV_ON && fid ~= -1
        fflush(fid);
        fclose(fid);
    end

    % plot dataset atual (limpa e desenha)
    cla(ax1); cla(ax2); cla(ax3);
    grid(ax1,'on'); grid(ax2,'on'); grid(ax3,'on');
    hold(ax1,'on'); hold(ax2,'on'); hold(ax3,'on');

    if ~isempty(all_t_hr), plot(ax1, all_t_hr, all_hr, 'LineWidth', 1.5); end
    if ~isempty(all_t_br), plot(ax2, all_t_br, all_br, 'LineWidth', 1.5); end

    if ~isempty(all_t_raw)
        plot(ax3, all_t_raw, all_hr_raw, 'LineWidth', 1.0);
        plot(ax3, all_t_raw, all_br_raw, 'LineWidth', 1.0);
        plot(ax3, all_t_raw, all_hr_f,   'LineWidth', 1.5);
        plot(ax3, all_t_raw, all_br_f,   'LineWidth', 1.5);
    end

    title(ax1, sprintf("%s | HR (tracker) | N=%d hop=%d | fs=%.2f", ds.name, CFG.N, CFG.hop, CFG.FS_TARGET));
    title(ax2, sprintf("%s | BR (tracker) | N=%d hop=%d | fs=%.2f", ds.name, CFG.N_BR, CFG.hop_BR, CFG.FS_TARGET));
    title(ax3, sprintf("%s | HR_COL=%d BR_COL=%d | HR BP=[%.3f %.3f]Hz | BR BP=[%.3f %.3f]Hz", ...
        ds.name, CFG.HR_COL, CFG.BR_COL, CFG.HR_FILT_WI, CFG.HR_FILT_WF, CFG.BR_FILT_WI, CFG.BR_FILT_WF));

    drawnow;
end

fprintf("\nDone.\n");

% ===================== LOCAL FUNCTIONS =====================

function tms = t_ms_from_tsec(t_sec, t0_ms)
    tms = double(t0_ms) + 1000.0 * double(t_sec);
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

function [sos, zi, err] = design_sos_cheby2(fs, order, wi, wf, rs)
    err = "";
    sos = [];
    zi  = [];

    fs    = double(fs);
    order = double(order);
    wi    = double(wi);
    wf    = double(wf);
    rs    = double(rs);

    if ~isfinite(fs) || fs <= 0
        err = "bad_fs"; return;
    end
    if mod(order,2) ~= 0 || order < 2
        err = "order_must_be_even"; return;
    end
    if ~isfinite(wi) || ~isfinite(wf) || wi <= 0 || wf <= wi
        err = "bad_wi_wf"; return;
    end
    if wf >= fs/2
        err = "wf>=nyq"; return;
    end

    try
        [z,p,k] = cheby2(order, rs, [wi wf]/(fs/2), 'bandpass');
        sos = zp2sos(z,p,k);
        zi  = sosfilt_zi(sos);
    catch ME
        err = "design_fail:" + string(ME.identifier);
        sos = []; zi = [];
    end
end

function trk = make_tracker(fs, N, hop, f_min_hz, prior_init_per_min, prior_min_per_min, prior_max_per_min, prior_sigma_hz, prior_alpha, ema_tau_sec, maxp)
    trk = struct();
    trk.fs = double(fs);
    trk.N = int32(N);
    trk.hop = int32(hop);
    trk.f_min_hz = double(f_min_hz);

    n = (0:double(N)-1)';
    trk.win = (0.21557895 ...
        - 0.41663158*cos(2*pi*n/double(N)) ...
        + 0.277263158*cos(4*pi*n/double(N)) ...
        - 0.083578947*cos(6*pi*n/double(N)) ...
        + 0.006947368*cos(8*pi*n/double(N)));

    nh = floor(double(N)/2) + 1;
    trk.f_bins = (0:nh-1)' * (trk.fs/double(N));   % rfft bins

    trk.k_min = max(2, int32(ceil(trk.f_min_hz * double(N) / trk.fs)) + 1); % +1 por DC (bin0) e garante >=2

    trk.prior_min = double(prior_min_per_min) / 60.0;
    trk.prior_max = double(prior_max_per_min) / 60.0;
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

    N = double(trk.N);

    % push no buffer (shift)
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

    if trk.len < trk.N
        return;
    end
    if trk.cnt < trk.hop
        return;
    end
    trk.cnt = int32(0);

    xw = trk.buf(:) .* trk.win(:);
    xw = xw - mean(xw);

    X = fft(xw, double(trk.N));
    nh = floor(double(trk.N)/2) + 1;
    P  = abs(X(1:nh)).^2;

    % faixa de busca (em bins 1-based de P)
    k_srch_min = max(double(trk.k_min), ceil(trk.prior_min * double(trk.N) / trk.fs) + 1);
    k_srch_max = min(numel(P), floor(trk.prior_max * double(trk.N) / trk.fs) + 1);

    if k_srch_max <= k_srch_min
        return;
    end

    % ===================== FIX: NADA de (P.*W)(a:b) =====================
    if isfinite(trk.prior_sigma_hz) && trk.prior_sigma_hz > 0
        W   = exp(-0.5 * ((trk.f_bins - trk.f_prev) / trk.prior_sigma_hz).^2);
        tmp = P(:) .* W(:);
        [~, kk] = max(tmp(k_srch_min:k_srch_max));
        k0 = (k_srch_min - 1) + kk;  % 1-based
    else
        [~, kk] = max(P(double(trk.k_min):end));
        k0 = double(trk.k_min) - 1 + kk; % 1-based
    end
    % ====================================================================

    kp    = parabolic_peak(P, k0);                 % fractional bin shift
    f_est = (double(k0-1) + kp) * trk.fs / double(trk.N);  % (k0-1) pq DC é bin 0

    if ~isfinite(f_est) || f_est <= 0
        return;
    end

    y_inst_per_min = 60.0 * f_est;

    if isnan(trk.ema_state)
        trk.ema_state = y_inst_per_min;
    else
        trk.ema_state = (1.0 - trk.ema_alpha)*trk.ema_state + trk.ema_alpha*y_inst_per_min;
    end

    f_clipped = min(max(f_est, trk.prior_min), trk.prior_max);
    trk.f_prev = (1.0 - trk.prior_alpha)*trk.f_prev + trk.prior_alpha*f_clipped;

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

function dp = parabolic_peak(P, k)
    dp = 0.0;
    if k <= 1 || k >= numel(P)
        return;
    end
    v   = log(P(k-1:k+1) + 1e-20);
    den = v(1) - 2*v(2) + v(3);
    if den == 0 || ~isfinite(den)
        return;
    end
    dp = 0.5 * (v(1) - v(3)) / den;
    dp = min(max(dp, -0.75), 0.75);
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
        x(:) = 0;
        return;
    end

    x(~isfinite(x)) = 0;
    if ~any(isfinite(x)), x(:) = 0; end
end