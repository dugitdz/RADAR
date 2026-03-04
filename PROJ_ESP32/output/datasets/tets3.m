clc; clear; close all;

%% ===================== PATHS / DATASETS =====================
BASE = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",    "radar", BASE+"phases.csv",      "polar", BASE+"POLARH2.txt")
    struct("name","H3",    "radar", BASE+"phases_raw.csv",  "polar", BASE+"POLARH3.txt")
    struct("name","H10",   "radar", BASE+"phase.csv",       "polar", BASE+"POLARH10.txt")
    struct("name","TESTE", "radar", BASE+"Radar_teste.csv", "polar", BASE+"Polar_teste.txt")
];

COL_T_MS  = 1;
COL_PHASE = 2;

%% ===================== CSV OUTPUT (igual Python) =====================
CSV_DIR         = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\HR_BR";
CSV_ON          = true;
CSV_FLUSH_EVERY = 50;

%% ===================== "PYTHON CFG" (mesma lógica) =====================
CFG = struct();

CFG.HR_COL = 1;   % 1=total 2=resp 3=heart
CFG.BR_COL = 2;

CFG.N      = 256;  % HR
CFG.hop    = 32;
CFG.N_BR   = 512;  % BR
CFG.hop_BR = 64;

CFG.FS_TARGET = 50.0;   % resample offline para fs fixo
CFG.GAP_THR   = 0.6;    % segundos
IMETH         = 'linear';

% HR bandpass causal (cheby2, ordem PAR)
CFG.HR_FILT_ORD = 4;
CFG.HR_FILT_WI  = 0.595;
CFG.HR_FILT_WF  = 3.05;
CFG.HR_FILT_RS  = 30.50;

% BR bandpass causal (cheby2, ordem PAR)
CFG.BR_FILT_ORD = 4;
CFG.BR_FILT_WI  = 0.10;
CFG.BR_FILT_WF  = 0.50;
CFG.BR_FILT_RS  = 30.0;

% prior HR
CFG.HR_PRIOR_INIT_BPM = 80.0;
CFG.HR_PRIOR_MIN_BPM  = 35.0;
CFG.HR_PRIOR_MAX_BPM  = 240.0;
CFG.HR_PRIOR_SIGMA_HZ = 0.40;
CFG.HR_PRIOR_ALPHA    = 0.05;
CFG.HR_EMA_TAU_SEC    = 4.0;

% prior BR
CFG.BR_PRIOR_INIT_BRPM = 15.0;
CFG.BR_PRIOR_MIN_BRPM  = 4.0;
CFG.BR_PRIOR_MAX_BRPM  = 40.0;
CFG.BR_PRIOR_SIGMA_HZ  = 0.10;
CFG.BR_PRIOR_ALPHA     = 0.10;
CFG.BR_EMA_TAU_SEC     = 6.0;

CFG.F_MIN_HZ = 0.05;

MAXP = 2000;

%% ===================== VALIDATION =====================
if CFG.N < 32 || CFG.hop < 1 || CFG.hop > CFG.N
    error("CFG inválido (HR): verifique N/hop");
end
if CFG.N_BR < 32 || CFG.hop_BR < 1 || CFG.hop_BR > CFG.N_BR
    error("CFG inválido (BR): verifique N_BR/hop_BR");
end
if ~ismember(CFG.HR_COL,[1 2 3]) || ~ismember(CFG.BR_COL,[1 2 3])
    error("CFG inválido: HR_COL/BR_COL devem ser 1,2,3");
end

%% ===================== DESIGN SOS FILTERS (MATLAB TOOLBOX) =====================
[sos_hr, zi_hr0] = design_sos_cheby2_bp(CFG.FS_TARGET, CFG.HR_FILT_ORD, CFG.HR_FILT_WI, CFG.HR_FILT_WF, CFG.HR_FILT_RS);
[sos_br, zi_br0] = design_sos_cheby2_bp(CFG.FS_TARGET, CFG.BR_FILT_ORD, CFG.BR_FILT_WI, CFG.BR_FILT_WF, CFG.BR_FILT_RS);

%% ===================== GUI (igual Python) =====================
fig = figure('Color','w','Name','Python-tracker -> MATLAB (dataset input)');
tl  = tiledlayout(fig,3,1,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl,1); grid(ax1,'on'); hold(ax1,'on'); ylabel(ax1,'HR (bpm)');
ax2 = nexttile(tl,2); grid(ax2,'on'); hold(ax2,'on'); ylabel(ax2,'BR (brpm)');
ax3 = nexttile(tl,3); grid(ax3,'on'); hold(ax3,'on'); ylabel(ax3,'signals'); xlabel(ax3,'t (s)');
legend(ax3, {'heart_raw','resp_raw','heart_filt','resp_filt'}, 'Location','best');

fprintf("\n=== Python-tracker traduzido p/ MATLAB (offline via arquivos) ===\n");
fprintf("FS_TARGET=%.2f | GAP_THR=%.2fs | HR N=%d hop=%d | BR N=%d hop=%d\n", ...
    CFG.FS_TARGET, CFG.GAP_THR, CFG.N, CFG.hop, CFG.N_BR, CFG.hop_BR);
fprintf("HR BP=[%.3f %.3f]Hz (cheby2 ord=%d Rs=%.1f)\n", CFG.HR_FILT_WI, CFG.HR_FILT_WF, CFG.HR_FILT_ORD, CFG.HR_FILT_RS);
fprintf("BR BP=[%.3f %.3f]Hz (cheby2 ord=%d Rs=%.1f)\n", CFG.BR_FILT_WI, CFG.BR_FILT_WF, CFG.BR_FILT_ORD, CFG.BR_FILT_RS);
fprintf("===============================================================\n");

%% ===================== MAIN LOOP (datasets) =====================
for d = 1:numel(DATASETS)
    ds = DATASETS(d);
    fprintf("\n---- %s ----\n", ds.name);

    % ---------- lê radar CSV ----------
    A = readmatrix(ds.radar);

    if size(A,2) >= 4
        t_ms  = A(:,1);
        total = A(:,2);
        resp  = A(:,3);
        heart = A(:,4);
    else
        t_ms  = A(:,COL_T_MS);
        ph    = A(:,COL_PHASE);
        total = ph; resp = ph; heart = ph;
        fprintf("Radar %d colunas -> usando phase como total/resp/heart.\n", size(A,2));
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
        fprintf("%s | radar vazio\n", ds.name);
        continue;
    end

    % normaliza tempo a partir do primeiro
    t0 = (t_ms - t_ms(1))/1000.0;
    [t0, iu] = unique(t0,'stable');
    t_ms  = t_ms(iu);
    total = total(iu);
    resp  = resp(iu);
    heart = heart(iu);

    % segmenta por gaps (no tempo original)
    segs = segment_by_gaps(t0, CFG.GAP_THR);

    % tracker structs (zerados por segmento)
    trk_hr = make_tracker(CFG.FS_TARGET, CFG.N,    CFG.hop,    CFG.F_MIN_HZ, ...
        CFG.HR_PRIOR_INIT_BPM, CFG.HR_PRIOR_MIN_BPM, CFG.HR_PRIOR_MAX_BPM, ...
        CFG.HR_PRIOR_SIGMA_HZ, CFG.HR_PRIOR_ALPHA, CFG.HR_EMA_TAU_SEC, MAXP);

    trk_br = make_tracker(CFG.FS_TARGET, CFG.N_BR, CFG.hop_BR, CFG.F_MIN_HZ, ...
        CFG.BR_PRIOR_INIT_BRPM, CFG.BR_PRIOR_MIN_BRPM, CFG.BR_PRIOR_MAX_BRPM, ...
        CFG.BR_PRIOR_SIGMA_HZ, CFG.BR_PRIOR_ALPHA, CFG.BR_EMA_TAU_SEC, MAXP);

    last_hr = NaN;
    last_br = NaN;

    % CSV por dataset
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

    % buffers para plot (dataset)
    t_raw    = [];
    heart_r  = [];
    breath_r = [];
    heart_f  = [];
    breath_f = [];

    t_hr = [];
    y_hr = [];
    t_br = [];
    y_br = [];

    dt = 1.0/CFG.FS_TARGET;

    % ---------- percorre segmentos ----------
    for s = 1:size(segs,1)
        i0 = segs(s,1); i1 = segs(s,2);

        ts = t0(i0:i1);
        tot_s = total(i0:i1);
        resp_s = resp(i0:i1);
        heart_s = heart(i0:i1);

        [ts, iu2] = unique(ts,'stable');
        tot_s   = tot_s(iu2);
        resp_s  = resp_s(iu2);
        heart_s = heart_s(iu2);

        if numel(ts) < 8
            continue;
        end

        % resample uniforme
        tnew = (ts(1):dt:ts(end))';
        tot_u   = interp1(ts, tot_s,   tnew, IMETH);
        resp_u  = interp1(ts, resp_s,  tnew, IMETH);
        heart_u = interp1(ts, heart_s, tnew, IMETH);

        ok2 = isfinite(tnew) & isfinite(tot_u) & isfinite(resp_u) & isfinite(heart_u);
        tnew    = tnew(ok2);
        tot_u   = force_finite_vector(tot_u(ok2));
        resp_u  = force_finite_vector(resp_u(ok2));
        heart_u = force_finite_vector(heart_u(ok2));

        if numel(tnew) < max(CFG.N_BR, CFG.N) + 2
            continue;
        end

        % reset "segment" igual Python: zi = zi0 * x_init, e reseta trackers
        x0_hr = pick_col(tot_u(1), resp_u(1), heart_u(1), CFG.HR_COL);
        x0_br = pick_col(tot_u(1), resp_u(1), heart_u(1), CFG.BR_COL);

        zi_hr = zi_hr0 * x0_hr;
        zi_br = zi_br0 * x0_br;

        trk_hr = tracker_reset(trk_hr);
        trk_br = tracker_reset(trk_br);

        % processa amostra-a-amostra (cópia do realtime)
        for i = 1:numel(tnew)
            t_sec = tnew(i);

            x_hr_raw = pick_col(tot_u(i), resp_u(i), heart_u(i), CFG.HR_COL);
            x_br_raw = pick_col(tot_u(i), resp_u(i), heart_u(i), CFG.BR_COL);

            [yhr, zi_hr] = sosfilt(sos_hr, x_hr_raw, zi_hr);
            [ybr, zi_br] = sosfilt(sos_br, x_br_raw, zi_br);
            yhr = double(yhr);
            ybr = double(ybr);

            [trk_hr, hr_est] = tracker_push(trk_hr, t_sec, yhr);
            [trk_br, br_est] = tracker_push(trk_br, t_sec, ybr);

            % buffers p/ plot
            t_raw(end+1,1)    = t_sec;
            heart_r(end+1,1)  = x_hr_raw;
            breath_r(end+1,1) = x_br_raw;
            heart_f(end+1,1)  = yhr;
            breath_f(end+1,1) = ybr;

            if ~isnan(hr_est)
                last_hr = hr_est;
                t_hr(end+1,1) = trk_hr.t_out(trk_hr.out_len);
                y_hr(end+1,1) = trk_hr.y_out(trk_hr.out_len);
            end
            if ~isnan(br_est)
                last_br = br_est;
                t_br(end+1,1) = trk_br.t_out(trk_br.out_len);
                y_br(end+1,1) = trk_br.y_out(trk_br.out_len);
            end

            % CSV: 1 linha por "pacote" (aqui: por amostra do resample)
            if CSV_ON
                t_ms_out = (t_ms(1) + 1000.0*t_sec);
                fprintf(fid, "%.3f,%.9g,%.9g,%.9g,%.9g,%.9g\n", ...
                    t_ms_out, tot_u(i), resp_u(i), heart_u(i), last_br, last_hr);
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

    % ---------- plot dataset ----------
    cla(ax1); cla(ax2); cla(ax3);
    grid(ax1,'on'); grid(ax2,'on'); grid(ax3,'on');
    hold(ax1,'on'); hold(ax2,'on'); hold(ax3,'on');

    if ~isempty(t_hr), plot(ax1, t_hr, y_hr, 'LineWidth', 1.5); end
    if ~isempty(t_br), plot(ax2, t_br, y_br, 'LineWidth', 1.5); end

    if ~isempty(t_raw)
        plot(ax3, t_raw, heart_r,  'LineWidth', 1.0);
        plot(ax3, t_raw, breath_r, 'LineWidth', 1.0);
        plot(ax3, t_raw, heart_f,  'LineWidth', 1.5);
        plot(ax3, t_raw, breath_f, 'LineWidth', 1.5);
    end

    title(ax1, sprintf("%s | HR tracker | N=%d hop=%d | fs=%.2f", ds.name, CFG.N, CFG.hop, CFG.FS_TARGET));
    title(ax2, sprintf("%s | BR tracker | N=%d hop=%d | fs=%.2f", ds.name, CFG.N_BR, CFG.hop_BR, CFG.FS_TARGET));
    title(ax3, sprintf("%s | HR_COL=%d BR_COL=%d | HR BP=[%.3f %.3f]Hz | BR BP=[%.3f %.3f]Hz", ...
        ds.name, CFG.HR_COL, CFG.BR_COL, CFG.HR_FILT_WI, CFG.HR_FILT_WF, CFG.BR_FILT_WI, CFG.BR_FILT_WF));

    drawnow;
end

fprintf("\nDone.\n");

%% ===================== FUNCTIONS =====================

function [sos, zi0] = design_sos_cheby2_bp(fs, ord_total, wi, wf, rs)
    fs = double(fs); wi = double(wi); wf = double(wf); rs = double(rs);
    ord_total = int32(ord_total);

    if mod(ord_total,2)~=0 || ord_total<2
        error("cheby2: ordem deve ser par (ord_total=%d)", ord_total);
    end
    if wf >= fs/2
        error("cheby2: wf>=Nyquist (wf=%.3f fs=%.3f)", wf, fs);
    end
    if wi <= 0 || wf <= wi
        error("cheby2: wi/wf inválidos (wi=%.3f wf=%.3f)", wi, wf);
    end

    % Python: cheby2(order//2, rs, ..., output="sos")
    % MATLAB: cheby2(N,Rs,Wn,'bandpass','sos') => N é ordem do protótipo (metade quando SOS de bandpass)
    Nproto = double(ord_total)/2;
    sos = cheby2(Nproto, rs, [wi wf]/(fs/2), 'bandpass', 'sos');
    zi0 = sosfilt_zi(sos);
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

function x = pick_col(total, resp, heart, col_idx)
    if col_idx == 1
        x = double(total);
    elseif col_idx == 2
        x = double(resp);
    else
        x = double(heart);
    end
end

function trk = make_tracker(fs, N, hop, f_min_hz, ...
    prior_init_per_min, prior_min_per_min, prior_max_per_min, ...
    prior_sigma_hz, prior_alpha, ema_tau_sec, maxp)

    trk = struct();
    trk.fs  = double(fs);
    trk.N   = int32(N);
    trk.hop = int32(hop);
    trk.f_min_hz = double(f_min_hz);

    % flattop (mesma do Python)
    n = (0:double(N)-1)';
    trk.win = (0.21557895 ...
        - 0.41663158*cos(2*pi*n/double(N)) ...
        + 0.277263158*cos(4*pi*n/double(N)) ...
        - 0.083578947*cos(6*pi*n/double(N)) ...
        + 0.006947368*cos(8*pi*n/double(N)));

    trk.f_bins = (0:floor(double(N)/2))' * (trk.fs/double(N)); % rfft freqs

    % k_min em bins (0-based no Python). No MATLAB vamos usar "k0bin" 0-based internamente.
    kmin0 = max(1, ceil(trk.f_min_hz * double(N) / trk.fs) - 1); % igual Python
    trk.k_min0 = double(kmin0);

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

    N = double(trk.N);

    % push (shift buffer)
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

    % window + demean
    xw = trk.buf(:) .* trk.win(:);
    xw = xw - mean(xw);

    % rFFT power
    X  = fft(xw, double(trk.N));
    nh = floor(double(trk.N)/2) + 1;
    P  = abs(X(1:nh)).^2;

    % search limits (bins 0-based -> converter)
    k_srch_min0 = max(trk.k_min0, ceil(trk.prior_min * double(trk.N) / trk.fs)); % 0-based
    k_srch_max0 = min(nh-1, floor(trk.prior_max * double(trk.N) / trk.fs));     % 0-based

    if k_srch_max0 <= k_srch_min0
        return;
    end

    % MATLAB indices 1-based
    i_min = int32(k_srch_min0 + 1);
    i_max = int32(k_srch_max0 + 1);

    if isfinite(trk.prior_sigma_hz) && trk.prior_sigma_hz > 0
        W = exp(-0.5 * ((trk.f_bins - trk.f_prev) / trk.prior_sigma_hz).^2);
        tmp = P(:) .* W(:);
        [~, kk] = max(tmp(double(i_min):double(i_max)));
        i0 = double(i_min) - 1 + kk; % 1-based index do MATLAB em P
    else
        [~, kk] = max(P(double(trk.k_min0+1):end));
        i0 = double(trk.k_min0+1) - 1 + kk;
    end

    % peak refine (parabolic, usando índice MATLAB 1-based)
    dp = parabolic_peak_log(P, i0);
    k0bin = (double(i0)-1) + dp;      % bin 0-based fracionário
    f_est = k0bin * trk.fs / double(trk.N);

    if ~isfinite(f_est) || f_est <= 0
        return;
    end

    y_inst = 60.0 * f_est;

    if isnan(trk.ema_state)
        trk.ema_state = y_inst;
    else
        trk.ema_state = (1.0 - trk.ema_alpha)*trk.ema_state + trk.ema_alpha*y_inst;
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

function dp = parabolic_peak_log(P, i)
    dp = 0.0;
    if i <= 1 || i >= numel(P)
        return;
    end
    v = log(P(i-1:i+1) + 1e-20);
    den = v(1) - 2*v(2) + v(3);
    if den == 0 || ~isfinite(den)
        return;
    end
    dp = 0.5 * (v(1) - v(3)) / den;
    dp = min(max(dp, -0.75), 0.75);
end