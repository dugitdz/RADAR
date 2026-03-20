%% Initialize
clearvars; close all; clc;

%% ===================== PATHS (4 DUPLAS) =====================
BASE = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",    "radar", BASE + "radar7.csv",   "polar", BASE + "polar7.txt")
    struct("name","H3",    "radar", BASE + "radar8.csv",   "polar", BASE + "polar8.txt")
    struct("name","H10",   "radar", BASE + "radar9.csv",   "polar", BASE + "polar9.txt")
    struct("name","TESTE", "radar", BASE + "radar10.csv",  "polar", BASE + "polar10.txt")
];

%% ===================== ALGO A (CSV HR DIRECT) =====================
ALGOA_MOVMEAN_N = 7;

%% ===================== ALGO B (PHASE -> DWT -> FFTPEAK) =====================
COL_T_MS  = 1;
COL_PHASE = 2;

gap_thr   = 0.5;
FS_TARGET = 25.0;
IMETH     = 'linear';
WRAP_ON   = 0;

DWT_WAVE  = 'db5';
DWT_LEVEL = 4;

winN = 64;
hopN = 16;
HOP_SEC = hopN/FS_TARGET;

FMIN_HZ = 0.9;
FMAX_HZ = 2.0;

CONF_MODE = "pmax_norm";
CONF_THR  = 0.015;

KEEPD = [4 3];

NFFT_BIG = 2048;
FINAL_MOVMEAN_SEC = 7;

%% ===================== PLOT HELPERS =====================
FILL_GAPS_FOR_PLOT = 1;
MAX_GAP_SEC_PLOT   = 4.0;

fprintf('\n================ RESULTADOS (RICARDO vs NOVO) ================\n');

for d = 1:numel(DATASETS)
    ds = DATASETS(d);

    radar_path = resolve_any_path(ds.radar);
    radar_label = get_file_stem(radar_path);

    polar_path = resolve_polar_path(ds.polar);
    [tP, hP] = read_txt_polar_flex(polar_path);
    if isempty(tP)
        fprintf('%s | Polar vazio. Pulando.\n', radar_label);
        continue;
    end

    [hA_onP, okA] = algoA_hr_direct_on_polar(radar_path, tP, ALGOA_MOVMEAN_N);
    hA_onP_plot = hA_onP;
    if okA && FILL_GAPS_FOR_PLOT
        hA_onP_plot = fill_gaps_limited(tP, hA_onP, MAX_GAP_SEC_PLOT);
    end

    [RMSEA, MAEA, CORRA, nA] = metrics_vs_polar(hA_onP, hP);

    [hB_onP, okB] = algoB_phase_on_polar(radar_path, tP, ...
        COL_T_MS, COL_PHASE, gap_thr, FS_TARGET, IMETH, WRAP_ON, ...
        DWT_WAVE, DWT_LEVEL, KEEPD, winN, hopN, FMIN_HZ, FMAX_HZ, ...
        CONF_MODE, CONF_THR, NFFT_BIG, FINAL_MOVMEAN_SEC);

    hB_onP_plot = hB_onP;
    if okB && FILL_GAPS_FOR_PLOT
        hB_onP_plot = fill_gaps_limited(tP, hB_onP, MAX_GAP_SEC_PLOT);
    end

    [RMSEB, MAEB, CORRB, nB] = metrics_vs_polar(hB_onP, hP);

    fprintf('%s | Ricardo: RMSE=%.3f MAE=%.3f CORR=%.3f N=%d | Novo: RMSE=%.3f MAE=%.3f CORR=%.3f N=%d\n', ...
        radar_label, RMSEA, MAEA, CORRA, nA, RMSEB, MAEB, CORRB, nB);

    figure('Color','w','Name', radar_label + " | HR (Ricardo,Novo) vs Polar");
    hold on; grid on;
    plot(tP, hP, 'r-', 'LineWidth', 1.6, 'DisplayName','POLAR');
    if okA
        plot(tP, hA_onP_plot, 'm-', 'LineWidth', 1.4, 'DisplayName', sprintf('Ricardo (movmean=%d)', ALGOA_MOVMEAN_N));
    end
    if okB
        plot(tP, hB_onP_plot, 'k-', 'LineWidth', 1.4, 'DisplayName', 'Novo (DWT+FFTpeak)');
    end
    xlabel('t (s)');
    ylabel('HR (bpm)');
    title(sprintf('%s | Ricardo: RMSE=%.2f MAE=%.2f Corr=%.2f | Novo: RMSE=%.2f MAE=%.2f Corr=%.2f', ...
        radar_label, RMSEA, MAEA, CORRA, RMSEB, MAEB, CORRB));
    legend('Location','best');
end

fprintf('==============================================================\n');

%% ===================== FUNCTIONS =====================

function p2 = resolve_any_path(p)
    p = string(p);
    if exist(p,'file') == 2
        p2 = char(p);
        return;
    end
    exts = [".csv" ".txt" ".CSV" ".TXT"];
    for e = exts
        if exist(p + e,'file') == 2
            p2 = char(p + e);
            return;
        end
    end
    p2 = char(p);
end

function p2 = resolve_polar_path(p)
    p2 = resolve_any_path(p);
end

function stem = get_file_stem(p)
    [~, stem, ~] = fileparts(char(p));
    stem = string(stem);
end

function [t_sec, HR] = read_txt_polar_flex(p)
    fid = fopen(p,'r');
    if fid == -1
        t_sec = [];
        HR = [];
        return;
    end
    raw = textscan(fid,'%s','Delimiter','\n');
    raw = raw{1};
    fclose(fid);

    if numel(raw) < 2
        t_sec = [];
        HR = [];
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
            HR(i)     = str2double(strtrim(parts(2)));
        end
    end

    try
        t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS");
    catch
        t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd HH:mm:ss.SSS");
    end

    t_sec = seconds(t_dt - t_dt(1));
    t_sec = t_sec(:);
    HR = HR(:);

    m = isfinite(t_sec) & isfinite(HR);
    t_sec = t_sec(m);
    HR = HR(m);

    [t_sec, iu] = unique(t_sec,'stable');
    HR = HR(iu);
end

function [h_onP, okA] = algoA_hr_direct_on_polar(radar_csv, tP, movN)
    okA = false;
    h_onP = nan(size(tP));

    try
        T = readtable(radar_csv);
    catch
        return;
    end

    if ~ismember("timestamp", string(T.Properties.VariableNames))
        return;
    end

    t_raw = T.timestamp;
    if isempty(t_raw) || numel(t_raw) < 5
        return;
    end

    if ismember("HR", string(T.Properties.VariableNames))
        y = T.HR;
    elseif ismember("heart", string(T.Properties.VariableNames))
        y = T.heart;
    else
        return;
    end

    dt_med = median(diff(t_raw), 'omitnan');
    if dt_med > 1
        t_sec = (t_raw - t_raw(1))/1000.0;
    else
        t_sec = (t_raw - t_raw(1));
    end

    m = isfinite(t_sec) & isfinite(y);
    t_sec = t_sec(m);
    y = y(m);

    if numel(t_sec) < 5
        return;
    end

    [t_sec, iu] = unique(t_sec,'stable');
    y = y(iu);

    y_sm = movmean(y, movN, 'omitnan', 'Endpoints','shrink');

    h_onP = interp1(t_sec, y_sm, tP, 'linear', NaN);
    okA = true;
end

function [RMSE, MAE, CORR, nValid] = metrics_vs_polar(h_onP, hP)
    m = isfinite(h_onP) & isfinite(hP);
    nValid = nnz(m);
    if nValid < 5
        RMSE = NaN;
        MAE = NaN;
        CORR = NaN;
        return;
    end
    e = h_onP(m) - hP(m);
    RMSE = sqrt(mean(e.^2,'omitnan'));
    MAE  = mean(abs(e),'omitnan');
    CORR = corr(h_onP(m), hP(m), 'Rows','complete');
end

function [h_onP, okB] = algoB_phase_on_polar(radar_path, tP, ...
    COL_T_MS, COL_PHASE, gap_thr, FS_TARGET, IMETH, WRAP_ON, ...
    DWT_WAVE, DWT_LEVEL, KEEPD, winN, hopN, FMIN_HZ, FMAX_HZ, ...
    CONF_MODE, CONF_THR, NFFT_BIG, FINAL_MOVMEAN_SEC)

    okB = false;
    h_onP = nan(size(tP));

    try
        A = readmatrix(radar_path);
    catch
        return;
    end
    if size(A,2) < max(COL_T_MS, COL_PHASE)
        return;
    end

    tpuro_ms = A(:, COL_T_MS);
    phase0   = A(:, COL_PHASE);

    tpuro_ms = tpuro_ms(:);
    phase0   = phase0(:);

    mabs = isfinite(tpuro_ms) & isfinite(phase0) & (abs(phase0) <= 10);
    tpuro_ms = tpuro_ms(mabs);
    phase0   = phase0(mabs);

    if numel(tpuro_ms) < 8
        return;
    end

    t0 = (tpuro_ms - tpuro_ms(1))/1000;
    [t0, ia] = unique(t0, 'stable');
    phase0   = phase0(ia);

    ph = force_finite_vector(double(phase0));
    ph = unwrap(ph);

    if WRAP_ON == 1
        x0 = wrap_phase(ph);
    else
        x0 = ph - mean(ph,'omitnan');
        x0 = force_finite_vector(x0);
    end

    segments = segment_by_gaps(t0, gap_thr);
    dt = 1/FS_TARGET;

    tseg_list  = {};
    hrseg_list = {};

    for s = 1:size(segments,1)
        i0 = segments(s,1);
        i1 = segments(s,2);
        ts = t0(i0:i1);
        xs = x0(i0:i1);

        [ts, iu] = unique(ts,'stable');
        xs = xs(iu);

        if numel(ts) < 8
            continue;
        end

        tnew = (ts(1):dt:ts(end))';
        xnew = interp1(ts, xs, tnew, IMETH);

        ok = isfinite(tnew) & isfinite(xnew);
        tnew = tnew(ok);
        xnew = xnew(ok);

        if numel(tnew) < winN
            continue;
        end

        xnew = force_finite_vector(xnew);

        xdiff = [0; diff(xnew)];
        xdiff = force_finite_vector(xdiff);
        xdiff = xdiff - mean(xdiff,'omitnan');

        [C,L] = wavedec(xdiff, DWT_LEVEL, DWT_WAVE);
        x_rec = idwt_keep_details(C, L, DWT_WAVE, DWT_LEVEL, KEEPD);
        x_hr  = sqrt_energy_match(x_rec, xdiff);

        [t_hr, hr_bpm, ~] = estimate_hr_by_fftpeak(tnew, x_hr, FS_TARGET, winN, hopN, ...
            FMIN_HZ, FMAX_HZ, CONF_MODE, CONF_THR, NFFT_BIG);

        if isempty(t_hr)
            continue;
        end

        if FINAL_MOVMEAN_SEC > 0
            HOP_SEC = hopN/FS_TARGET;
            W = max(1, round(FINAL_MOVMEAN_SEC / HOP_SEC));
            hr_bpm = movmean(hr_bpm, W, 'omitnan', 'Endpoints','shrink');
        end

        tseg_list{end+1}  = t_hr(:);
        hrseg_list{end+1} = hr_bpm(:);
    end

    [tR, hR] = concat_segments(tseg_list, hrseg_list);
    if isempty(tR)
        return;
    end

    m = isfinite(tR) & isfinite(hR);
    if nnz(m) < 5
        return;
    end

    h_onP = interp1(tR(m), hR(m), tP, 'linear', NaN);
    okB = true;
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

function [t_all, h_all] = concat_segments(tseg_list, hseg_list)
    t_all = [];
    h_all = [];
    for k = 1:numel(tseg_list)
        tS = tseg_list{k}(:);
        hS = hseg_list{k}(:);

        m = isfinite(tS) & isfinite(hS);
        tS = tS(m);
        hS = hS(m);

        if numel(tS) < 2
            continue;
        end

        [tS, iu] = unique(tS, 'stable');
        hS = hS(iu);

        t_all = [t_all; tS];
        h_all = [h_all; hS];
    end

    if isempty(t_all)
        return;
    end

    [t_all, iu] = unique(t_all, 'stable');
    h_all = h_all(iu);
end

function x = wrap_phase(ph)
    ph = double(ph(:));
    ph = force_finite_vector(ph);
    x = mod(ph + pi, 2*pi) - pi;
    x = x - mean(x,'omitnan');
end

function x = force_finite_vector(x)
    x = double(x(:));
    ok = isfinite(x);
    if all(ok)
        return;
    end
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
    if ~any(isfinite(x))
        x(:) = 0;
    end
end

function x_rec = idwt_keep_details(C, L, wname, Nlevel, keepD)
    keepD = unique(keepD(:))';
    keepD = keepD(keepD>=1 & keepD<=Nlevel);

    Cnew = zeros(size(C));
    idx = 1;

    lenA = L(1);
    idx = idx + lenA;

    for j = 1:Nlevel
        lev  = Nlevel - j + 1;
        lenD = L(j+1);
        i1 = idx;
        i2 = idx + lenD - 1;

        if any(keepD == lev)
            Cnew(i1:i2) = C(i1:i2);
        end

        idx = i2 + 1;
    end

    x_rec = waverec(Cnew, L, wname);
    x_rec = x_rec(:);
    x_rec = x_rec - mean(x_rec,'omitnan');
    x_rec = fillmissing(x_rec,'linear','EndValues','nearest');
    x_rec(~isfinite(x_rec)) = 0;
end

function [t_ctr, hr_bpm, conf_vec] = estimate_hr_by_fftpeak(t, x, Fs, winN, hopN, fmin_hz, fmax_hz, conf_mode, conf_thr, nfft_big)
    t = t(:);
    x = double(x(:));
    x = x - mean(x,'omitnan');
    x = fillmissing(x,'linear','EndValues','nearest');
    x(~isfinite(x)) = 0;

    N = numel(x);
    if N < winN
        t_ctr = [];
        hr_bpm = [];
        conf_vec = [];
        return;
    end

    idx0 = 1:hopN:(N-winN+1);
    nt = numel(idx0);

    t_ctr    = zeros(nt,1);
    hr_bpm   = nan(nt,1);
    conf_vec = nan(nt,1);

    for i = 1:nt
        i1 = idx0(i);
        seg = x(i1:i1+winN-1);
        seg = seg - mean(seg,'omitnan');

        nfft = max(nfft_big, 2^nextpow2(winN));
        X = fft(seg, nfft);
        P = abs(X).^2;

        f = (0:nfft-1)'*(Fs/nfft);
        nh = floor(nfft/2)+1;
        P = P(1:nh);
        f = f(1:nh);

        m = (f >= fmin_hz) & (f <= min(fmax_hz, Fs/2));
        t_ctr(i) = t(i1 + floor(winN/2));

        if nnz(m) < 5
            continue;
        end

        Pb = P(m);
        if ~any(isfinite(Pb)) || max(Pb) <= 0
            continue;
        end

        [pmax, kmax] = max(Pb);
        fb = f(m);

        if conf_mode == "pmax_norm"
            conf = pmax / (sum(Pb) + eps);
        else
            conf = 1;
        end
        conf_vec(i) = conf;

        if conf >= conf_thr
            hr_bpm(i) = 60*fb(kmax);
        end
    end
end

function y = sqrt_energy_match(x_rec, x_ref)
    x_rec = double(x_rec(:));
    x_ref = double(x_ref(:));
    Er = sum(x_rec.^2, 'omitnan');
    E0 = sum(x_ref.^2, 'omitnan');
    if ~isfinite(Er) || Er <= 0 || ~isfinite(E0) || E0 <= 0
        y = x_rec;
        return;
    end
    g = sqrt(E0 / Er);
    y = g * x_rec;
end

function yplot = fill_gaps_limited(t, y, max_gap_sec)
    t = t(:);
    y = y(:);
    yplot = y;

    if numel(t) < 3
        return;
    end

    isn = ~isfinite(yplot);
    if ~any(isn)
        return;
    end

    d = diff([0; isn; 0]);
    iStart = find(d == 1);
    iEnd   = find(d == -1) - 1;

    for k = 1:numel(iStart)
        a = iStart(k);
        b = iEnd(k);

        if a == 1 || b == numel(yplot)
            continue;
        end

        gap_sec = t(b) - t(a);
        if gap_sec <= max_gap_sec
            yplot(a:b) = interp1([t(a-1); t(b+1)], [yplot(a-1); yplot(b+1)], t(a:b), 'linear');
        end
    end
end