clc; clear; close all;

BASE = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",  "radar", BASE + "phases.csv",     "polar", BASE + "POLARH2.txt")
    struct("name","H3",  "radar", BASE + "phases_raw.csv", "polar", BASE + "POLARH3.txt")
    struct("name","H10", "radar", BASE + "phase.csv",      "polar", BASE + "POLARH10.txt")
];

COL_T_MS  = 1;
COL_PHASE = 2;

gap_thr   = 0.5;
FS_TARGET = 20.0;
IMETH     = 'linear';
WRAP_ON   = 0;

DWT_WAVE  = 'db5';
DWT_LEVEL = 4;

HR_LEAF   = "D3";        % "D3" (paper) ou "D2D3" (teste)

WIN_SEC   = 128/FS_TARGET;
HOP_SEC   = 32/FS_TARGET;
USE_WELCH = 1;

FMIN_HZ   = 0.8;
FMAX_HZ   = 2.0;

CONF_MODE = "pmax_norm"; % "pmax_norm" ou "peak_over_median"
CONF_THR  = 0.01;

FINAL_MOVMEAN_SEC = 7;

SPEC_BPM_MAX = 300;
SPEC_SMOOTH_BINS = 1;

fig = figure('Color','w','Name','HR vs Polar + Espectro');
tg  = uitabgroup(fig);

nD = numel(DATASETS);

fprintf('\n================ RESULTADOS (RMSE / MAE / CORR) ================\n');

for d = 1:nD
    ds = DATASETS(d);

    A = readmatrix(ds.radar);
    tpuro_ms = A(:, COL_T_MS);
    phase0   = A(:, COL_PHASE);

    tpuro_ms = tpuro_ms(:);
    phase0   = phase0(:);

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

    [tP, hP] = read_txt_polar_flex(ds.polar);
    tP = tP(:); hP = hP(:);

    if isempty(tP)
        fprintf('%s | Polar vazio. Pulando.\n', ds.name);
        continue;
    end

    dt = 1/FS_TARGET;
    tseg_list  = {};
    hrseg_list = {};

    best_len  = -inf;
    best_sum  = [];
    best_tsum = [];

    for s = 1:size(segments,1)
        i0 = segments(s,1); i1 = segments(s,2);
        ts = t0(i0:i1);
        xs = x0(i0:i1);

        [ts, iu] = unique(ts, 'stable');
        xs = xs(iu);

        if numel(ts) < 8
            continue;
        end

        tnew = (ts(1):dt:ts(end))';
        xnew = interp1(ts, xs, tnew, IMETH);

        ok = isfinite(tnew) & isfinite(xnew);
        tnew = tnew(ok);
        xnew = xnew(ok);

        if numel(tnew) < 32
            continue;
        end

        xnew = force_finite_vector(xnew);

        xdiff = [0; diff(xnew)];
        xdiff = force_finite_vector(xdiff);
        xdiff = xdiff - mean(xdiff,'omitnan');

        [C,L] = wavedec(xdiff, DWT_LEVEL, DWT_WAVE);

        D = cell(DWT_LEVEL,1);
        for lev = 1:DWT_LEVEL
            D{lev} = force_finite_vector(wrcoef('d', C, L, DWT_WAVE, lev));
        end
        A4 = force_finite_vector(wrcoef('a', C, L, DWT_WAVE, DWT_LEVEL));

        if HR_LEAF == "D2D3"
            x_leaf = force_finite_vector(D{4} + D{3});
        else
            x_leaf = force_finite_vector(D{3});
        end

        x_hr = x_leaf;
        x_hr = sqrt_energy_match(x_hr, x_leaf);

        seg_len = tnew(end) - tnew(1);
        if seg_len > best_len
            best_len  = seg_len;
            best_sum  = x_hr(:);
            best_tsum = tnew(:);
        end

        [t_hr, hr_bpm] = estimate_hr_by_psd(tnew, x_hr, FS_TARGET, WIN_SEC, HOP_SEC, USE_WELCH, FMIN_HZ, FMAX_HZ, CONF_MODE, CONF_THR);

        if isempty(t_hr)
            continue;
        end

        if FINAL_MOVMEAN_SEC > 0
            W = max(1, round(FINAL_MOVMEAN_SEC / HOP_SEC));
            hr_bpm = movmean(hr_bpm, W, 'Endpoints','shrink');
        end

        tseg_list{end+1}  = t_hr(:);
        hrseg_list{end+1} = hr_bpm(:);
    end

    [tR, hR] = concat_segments(tseg_list, hrseg_list);
    if isempty(tR)
        fprintf('%s | Sem HR estimado.\n', ds.name);
        continue;
    end

    hR_onP = interp1(tR, hR, tP, 'linear', NaN);

    mE = isfinite(hR_onP) & isfinite(hP);
    nValid = nnz(mE);

    if nValid >= 5
        e = hR_onP(mE) - hP(mE);
        RMSE = sqrt(mean(e.^2,'omitnan'));
        MAE  = mean(abs(e),'omitnan');
        CORR = corr(hR_onP(mE), hP(mE), 'Rows','complete');
    else
        RMSE = NaN; MAE = NaN; CORR = NaN;
    end

    fprintf('%s | RMSE=%.4f | MAE=%.4f | CORR=%.4f | N=%d\n', ds.name, RMSE, MAE, CORR, nValid);

    tab = uitab(tg, 'Title', ds.name);
    tl = tiledlayout(tab, 2, 1, 'Padding','compact', 'TileSpacing','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    plot(ax1, tP, hP,     'c--', 'LineWidth', 1.6);
    plot(ax1, tP, hR_onP, 'k-',  'LineWidth', 1.4);
    xlabel(ax1,'t (s)'); ylabel(ax1,'HR (bpm)');
    title(ax1, sprintf('%s | %s | Fs=%.1fHz | WIN=%.1fs HOP=%.1fs | THR=%.3f | RMSE=%.3f MAE=%.3f Corr=%.3f | N=%d', ...
        ds.name, HR_LEAF, FS_TARGET, WIN_SEC, HOP_SEC, CONF_THR, RMSE, MAE, CORR, nValid));
    legend(ax1,'POLAR','RADAR (interp@POLAR)','Location','best');

    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    if ~isempty(best_sum)
        [bpm, Pdb_sum] = spectrum_db(best_sum, FS_TARGET, SPEC_SMOOTH_BINS);
        plot(ax2, bpm, Pdb_sum, 'LineWidth', 1.3);
        xlim(ax2,[0 SPEC_BPM_MAX]);
        ylim(ax2,[-70 5]);
        xlabel(ax2,'BPM');
        ylabel(ax2,'Potência (dB rel. ao pico)');
        title(ax2, sprintf('Espectro (FFT) do HR no maior segmento (%.1fs)', best_len));
    else
        text(ax2, 0.1, 0.5, 'Sem segmento válido para espectro.', 'Units','normalized');
        axis(ax2,'off');
    end

    t_all = t0(:);
    x_all = x0(:);
    [t_all, iu] = unique(t_all, 'stable');
    x_all = x_all(iu);

    tnew_all = (t_all(1):dt:t_all(end))';
    xnew_all = interp1(t_all, x_all, tnew_all, IMETH);

    ok = isfinite(tnew_all) & isfinite(xnew_all);
    tnew_all = tnew_all(ok);
    xnew_all = xnew_all(ok);
    xnew_all = force_finite_vector(xnew_all);

    xdiff_all = [0; diff(xnew_all)];
    xdiff_all = force_finite_vector(xdiff_all);
    xdiff_all = xdiff_all - mean(xdiff_all,'omitnan');

    [C_all, L_all] = wavedec(xdiff_all, DWT_LEVEL, DWT_WAVE);

    D_all = cell(DWT_LEVEL,1);
    for lev = 1:DWT_LEVEL
        D_all{lev} = force_finite_vector(wrcoef('d', C_all, L_all, DWT_WAVE, lev));
    end
    A_all = force_finite_vector(wrcoef('a', C_all, L_all, DWT_WAVE, DWT_LEVEL));

    if HR_LEAF == "D2D3"
        leaf_all = force_finite_vector(D_all{2} + D_all{3});
    else
        leaf_all = force_finite_vector(D_all{3});
    end

    SUM_all = akf_smooth_adaptive(leaf_all);
    SUM_all = sqrt_energy_match(SUM_all, leaf_all);
    SUM_all = force_finite_vector(SUM_all);

    nrows = DWT_LEVEL + 2;
    fLeaves = figure('Color','w', 'Name', sprintf('%s - FFT das folhas (%s L%d)', ds.name, DWT_WAVE, DWT_LEVEL));
    tlv = tiledlayout(fLeaves, nrows, 1, 'Padding','compact', 'TileSpacing','compact');

    ax = nexttile(tlv,1); hold(ax,'on'); grid(ax,'on');
    [bpmA, PdbA] = spectrum_db(A_all, FS_TARGET, SPEC_SMOOTH_BINS);
    plot(ax, bpmA, PdbA, 'LineWidth', 1.2);
    xlim(ax,[0 SPEC_BPM_MAX]);
    ylim(ax,[-70 5]);
    ylabel(ax, sprintf('FFT A%d', DWT_LEVEL));
    title(ax, sprintf('%s | FFT (dB rel. pico) | HR=%s', ds.name, HR_LEAF));

    rr = 2;
    for lev = DWT_LEVEL:-1:1
        ax = nexttile(tlv,rr); hold(ax,'on'); grid(ax,'on');
        [bpmD, PdbD] = spectrum_db(D_all{lev}, FS_TARGET, SPEC_SMOOTH_BINS);
        plot(ax, bpmD, PdbD, 'LineWidth', 1.0);
        xlim(ax,[0 SPEC_BPM_MAX]);
        ylim(ax,[-70 5]);
        ylabel(ax, sprintf('FFT D%d', lev));
        rr = rr + 1;
    end

    ax = nexttile(tlv,nrows); hold(ax,'on'); grid(ax,'on');
    [bpmS, PdbS] = spectrum_db(SUM_all, FS_TARGET, SPEC_SMOOTH_BINS);
    plot(ax, bpmS, PdbS, 'LineWidth', 1.3);
    xlim(ax,[0 SPEC_BPM_MAX]);
    ylim(ax,[-70 5]);
    ylabel(ax,'FFT HR (suavizado)');
    xlabel(ax,'BPM');
end

fprintf('==============================================================\n\n');

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
        tS = tseg_list{k};
        hS = hseg_list{k};
        if numel(tS) < 2
            continue;
        end
        [tS, iu] = unique(tS, 'stable');
        hS = hS(iu);
        t_all = [t_all; tS];
        h_all = [h_all; hS];
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

function [t_ctr, hr_bpm] = estimate_hr_by_psd(t, x, Fs, win_sec, hop_sec, useWelch, fmin_hz, fmax_hz, conf_mode, conf_thr)
    t = t(:); x = double(x(:));
    x = x - mean(x,'omitnan');
    x = fillmissing(x,'linear','EndValues','nearest');
    x(~isfinite(x)) = 0;

    N = numel(x);
    winN = max(16, round(win_sec * Fs));
    hopN = max(1,  round(hop_sec * Fs));

    if N < winN
        t_ctr = [];
        hr_bpm = [];
        return;
    end

    idx0 = 1:hopN:(N-winN+1);
    nt = numel(idx0);

    t_ctr  = zeros(nt,1);
    hr_bpm = nan(nt,1);

    for i = 1:nt
        i1 = idx0(i);
        seg = x(i1:i1+winN-1);
        seg = seg - mean(seg,'omitnan');

        if useWelch == 1
            nfft  = max(256, 2^nextpow2(winN));
            nover = round(0.5*winN);
            [P,f] = pwelch(seg, hamming(winN,'periodic'), nover, nfft, Fs);
            P = P(:); f = f(:);
        else
            nfft = 2^nextpow2(winN);
            X = fft(seg, nfft);
            P = abs(X).^2;
            f = (0:nfft-1)'*(Fs/nfft);
            nh = floor(nfft/2)+1;
            P = P(1:nh);
            f = f(1:nh);
        end

        m = (f >= fmin_hz) & (f <= min(fmax_hz, Fs/2));
        t_ctr(i) = t(i1 + floor(winN/2));

        if nnz(m) < 5
            continue;
        end

        Pb = P(m);
        fb = f(m);

        if ~any(isfinite(Pb)) || max(Pb) <= 0
            continue;
        end

        [pmax, kmax] = max(Pb);

        ok = false;

        if conf_mode == "pmax_norm"
            denom = sum(Pb) + eps;
            conf = pmax / denom;
            ok = (conf >= conf_thr);

        elseif conf_mode == "peak_over_median"
            med = median(Pb(Pb>0));
            if ~isfinite(med) || med <= 0, med = eps; end
            conf = pmax / med;
            ok = (conf >= conf_thr);

        else
            ok = true;
        end

        if ok
            hr_bpm(i) = 60*fb(kmax);
        end
    end
end

function [bpm, Pdb] = spectrum_db(x, Fs, smooth_bins)
    x = double(x(:));
    x = x - mean(x,'omitnan');
    x = fillmissing(x,'linear','EndValues','nearest');
    x(~isfinite(x)) = 0;

    N = numel(x);
    nfft = 2^nextpow2(N);
    X = fft(x, nfft);
    P = abs(X).^2;

    f = (0:nfft-1)'*(Fs/nfft);
    nh = floor(nfft/2)+1;
    f = f(1:nh);
    P = P(1:nh);

    bpm = 60*f;

    Pref = max(P) + eps;
    Pdb = 10*log10(P/Pref + eps);

    if smooth_bins > 1
        Pdb = movmean(Pdb, smooth_bins, 'omitnan');
    end
end

function y = akf_smooth_adaptive(x)
    x = double(x(:));
    x = x - mean(x,'omitnan');
    x = fillmissing(x,'linear','EndValues','nearest');
    x(~isfinite(x)) = 0;

    N = numel(x);
    if N < 3
        y = x;
        return;
    end

    A = 1; H = 1;

    q = var(diff(x), 'omitnan');
    if ~isfinite(q) || q <= 0
        q = 1e-6;
    end
    Q = q;

    winR = min(40, max(10, floor(N/10)));

    xhat = 0;
    P = 1;

    y = zeros(N,1);

    for k = 1:N
        z = x(k);

        xhatp = A * xhat;
        Pp    = A * P * A + Q;

        i0 = max(1, k-winR+1);
        rr = var(x(i0:k) - mean(x(i0:k),'omitnan'), 'omitnan');
        if ~isfinite(rr) || rr <= 0
            rr = 1e-4;
        end
        R = rr;

        K  = Pp * H / (H * Pp * H + R);

        xhat = xhatp + K * (z - H * xhatp);
        P    = (1 - K * H) * Pp;

        y(k) = xhat;
    end

    y = y - mean(y,'omitnan');
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