clc; clear; close all;

BASE = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",    "radar", BASE + "phases.csv",        "polar", BASE + "POLARH2.txt")
    struct("name","H3",    "radar", BASE + "phases_raw.csv",    "polar", BASE + "POLARH3.txt")
    struct("name","H10",   "radar", BASE + "phase.csv",         "polar", BASE + "POLARH10.txt")
    struct("name","TESTE", "radar", BASE + "Radar_teste.csv",   "polar", BASE + "Polar_teste.txt")
];

COL_T_MS  = 1;
COL_PHASE = 2;

gap_thr   = 0.5;
FS_TARGET = 25.0;
IMETH     = 'linear';
WRAP_ON   = 0;

DWT_WAVE  = 'db5';
DWT_LEVEL = 4;

% ======= CONFIG ESCOLHIDA (a sua linha) =======
winN = 64;
hopN = 16;
WIN_SEC = winN/FS_TARGET;
HOP_SEC = hopN/FS_TARGET;

FMIN_HZ = 0.9;
FMAX_HZ = 2.0;

CONF_MODE = "pmax_norm";   % conf = pmax/sum(Pb)
CONF_THR  = 0.006;

KEEPD = [4 3];             % "[4 3]"

NFFT_BIG = 4096;
FINAL_MOVMEAN_SEC = 7;

% ===== plots =====
SPEC_BPM_MAX = 300;
SPEC_SMOOTH_BINS = 1;

CONF_DEBUG = 1;

FILL_GAPS_FOR_PLOT = 1;
MAX_GAP_SEC_PLOT   = 4.0;

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

    conf_all = [];
    conf_ok  = 0;
    conf_tot = 0;

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

        if numel(tnew) < winN
            continue;
        end

        xnew = force_finite_vector(xnew);

        % phase difference
        xdiff = [0; diff(xnew)];
        xdiff = force_finite_vector(xdiff);
        xdiff = xdiff - mean(xdiff,'omitnan');

        % DWT
        [C,L] = wavedec(xdiff, DWT_LEVEL, DWT_WAVE);

        % IDWT correta por seleção de detalhes (KEEPD)
        x_rec = idwt_keep_details(C, L, DWT_WAVE, DWT_LEVEL, KEEPD);

        % square-root normalization (energia)
        x_hr  = sqrt_energy_match(x_rec, xdiff);

        seg_len = tnew(end) - tnew(1);
        if seg_len > best_len
            best_len  = seg_len;
            best_sum  = x_hr(:);
            best_tsum = tnew(:);
        end

        [t_hr, hr_bpm, conf_vec] = estimate_hr_by_fftpeak(tnew, x_hr, FS_TARGET, winN, hopN, ...
                                                         FMIN_HZ, FMAX_HZ, CONF_MODE, CONF_THR, NFFT_BIG);

        if ~isempty(conf_vec)
            conf_all = [conf_all; conf_vec(:)];
            conf_tot = conf_tot + numel(conf_vec);
            conf_ok  = conf_ok  + nnz(isfinite(hr_bpm));
        end

        if isempty(t_hr)
            continue;
        end

        if FINAL_MOVMEAN_SEC > 0
            W = max(1, round(FINAL_MOVMEAN_SEC / HOP_SEC));
            hr_bpm = movmean(hr_bpm, W, 'omitnan', 'Endpoints','shrink');
        end

        tseg_list{end+1}  = t_hr(:);
        hrseg_list{end+1} = hr_bpm(:);
    end

    if CONF_DEBUG == 1
        if isempty(conf_all)
            fprintf('%s | CONF: sem janelas válidas.\n', ds.name);
        else
            c = conf_all(isfinite(conf_all));
            if isempty(c)
                fprintf('%s | CONF: tudo NaN.\n', ds.name);
            else
                fprintf('%s | CONF(pmax/sum): min=%.4g med=%.4g mean=%.4g max=%.4g | ok=%d/%d (thr=%.4g)\n', ...
                    ds.name, min(c), median(c), mean(c), max(c), conf_ok, conf_tot, CONF_THR);
            end
        end
    end

    [tR, hR] = concat_segments(tseg_list, hrseg_list);
    if isempty(tR)
        fprintf('%s | Sem HR estimado.\n', ds.name);
        continue;
    end

    % interp robusto (remove NaNs)
    m = isfinite(tR) & isfinite(hR);
    hR_onP = interp1(tR(m), hR(m), tP, 'linear', NaN);

    % apenas para plot (preenche buracos pequenos)
    hR_onP_plot = hR_onP;
    if FILL_GAPS_FOR_PLOT == 1
        hR_onP_plot = fill_gaps_limited(tP, hR_onP, MAX_GAP_SEC_PLOT);
    end

    % métricas com série crua (sem preencher)
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
    plot(ax1, tP, hP,          'c--', 'LineWidth', 1.6);
    plot(ax1, tP, hR_onP_plot, 'k-',  'LineWidth', 1.4);
    xlabel(ax1,'t (s)'); ylabel(ax1,'HR (bpm)');
    title(ax1, sprintf('%s | KEEPD=%s | Fs=%.1fHz | win=%d hop=%d | f=[%.1f..%.1f]Hz | thr=%.4g | RMSE=%.3f MAE=%.3f Corr=%.3f | N=%d', ...
        ds.name, mat2str(KEEPD), FS_TARGET, winN, hopN, FMIN_HZ, FMAX_HZ, CONF_THR, RMSE, MAE, CORR, nValid));
    legend(ax1,'POLAR','RADAR (interp@POLAR, gap-fill p/ plot)','Location','best');

    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    if ~isempty(best_sum)
        [bpm, Pn_sum] = spectrum_lin(best_sum, FS_TARGET, SPEC_SMOOTH_BINS);
        plot(ax2, bpm, Pn_sum, 'LineWidth', 1.3);
        xlim(ax2,[0 SPEC_BPM_MAX]);
        ylim(ax2,[0 1.05]);
        xlabel(ax2,'BPM');
        ylabel(ax2,'Potência (rel. ao pico)');
        title(ax2, sprintf('Espectro (FFT linear) do HR no maior segmento (%.1fs)', best_len));
    else
        text(ax2, 0.1, 0.5, 'Sem segmento válido para espectro.', 'Units','normalized');
        axis(ax2,'off');
    end
end

fprintf('==============================================================\n\n');

% ===================== FUNCTIONS =====================

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

    lenA = L(1);   % A_N
    idx = idx + lenA;

    % D_N ... D_1
    for j = 1:Nlevel
        lev = Nlevel - j + 1;
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
    t = t(:); x = double(x(:));
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

    t_ctr   = zeros(nt,1);
    hr_bpm  = nan(nt,1);
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

function [bpm, Pn] = spectrum_lin(x, Fs, smooth_bins)
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
    Pn = P / Pref;

    if smooth_bins > 1
        Pn = movmean(Pn, smooth_bins, 'omitnan');
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

function yplot = fill_gaps_limited(t, y, max_gap_sec)
    t = t(:); y = y(:);
    yplot = y;

    if numel(t) < 3
        return;
    end

    isn = ~isfinite(yplot);
    if ~any(isn)
        return;
    end

    d = diff([0; isn; 0]);
    iStart = find(d== 1);
    iEnd   = find(d==-1) - 1;

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