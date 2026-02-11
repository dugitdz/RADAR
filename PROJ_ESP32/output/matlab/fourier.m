clc; clear; close all;

%% ========================== PATHS (3 DUPLAS) ==========================
BASE = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",  "radar", BASE+"phases.csv",      "polar", BASE+"POLARH2.txt")
    struct("name","H3",  "radar", BASE+"phases_raw.csv",  "polar", BASE+"POLARH3.txt")
    struct("name","H10", "radar", BASE+"phase.csv",       "polar", BASE+"POLARH10.txt")
];

%% ===================== FIXOS =====================
GAP_THR_SEC = 0.6;      % segmentação por buracos (s)
FS_TARGET   = 25.0;     % Hz (resample uniforme)

COL_T_MS  = 1;          % ms
COL_HEART = 2;          % coluna do sinal (fase/heart)

DT_GRID       = 0.1;    % grid comum para métricas
SMOOTH_HR_SEC = 7.0;    % suavização do HR no grid comum
F_MIN_HZ      = 0.05;   % ignora DC/baixíssima freq no pico

MIN_POINTS_AFTER_RESAMPLE = 32;
MIN_METRIC_SAMPLES        = 10;

%% ===================== CONFIG (FFT) =====================
CFG.N    = 128;
CFG.hop  = 32;
CFG.WRAP = 0;

CFG.win   = "flattop";   % FIXO
CFG.ftype = "cheby2";    % FIXO
CFG.ford  = 4;           % ordem final (par) -> cheby2 usa n=ford/2
CFG.wi    = 0.60;        % Hz
CFG.wf    = 3.0;        % Hz

CFG.harm_on    = 1;
CFG.harm_max   = 3;
CFG.harm_tol   = 1;
CFG.harm_ratio = 0.30;

%% ===================== GATING (BP_THEN_GATE) =====================
GATE_ON     = 1;

HR_MIN_BPM  = 55;
HR_MAX_BPM  = 210;

DELTA_DB    = 0;        % capPow = Ppk_in * 10^(dB/10)
SMOOTH_BINS = 3;        % suaviza P(k) só pra estimar Ppk_in (0/1 desliga)
GATE_FLOOR  = 0.0;      % piso absoluto de potência (0 desliga)
GATE_MODE   = "CAP";    % "CAP" | "ZERO"

% Debug: reconstrução BP+GATE em tempo (overlap-add)
DEBUG_RECON_ON = true;

% Debug do gate: plota espectro 1 frame (pra ver máscara e teto)
DEBUG_GATE_PLOT = true;
DEBUG_GATE_FRAME = 10;      % frame (1..Nframes)
DEBUG_GATE_DATASET = "H2";  % H2/H3/H10

fprintf('\n==================== CONFIG DEFINITIVA ====================\n');
fprintf('N=%d hop=%d WRAP=%d win=%s | filt=%s ord=%d wi=%.2f wf=%.2f | harm=%d (max=%d tol=%d ratio=%.2f)\n', ...
    CFG.N, CFG.hop, CFG.WRAP, CFG.win, CFG.ftype, CFG.ford, CFG.wi, CFG.wf, ...
    CFG.harm_on, CFG.harm_max, CFG.harm_tol, CFG.harm_ratio);

fprintf('GATE_ON=%d | HR=[%d..%d] bpm | DELTA_DB=%g dB | SMOOTH_BINS=%d | GATE_FLOOR=%g | MODE=%s\n', ...
    GATE_ON, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR, GATE_MODE);

%% ===================== PREP (FILTRO + JANELA) =====================
% janela do FFT tracking
if CFG.win == "flattop"
    WIN = flattopwin(CFG.N,'periodic');
elseif CFG.win == "hann"
    WIN = hann(CFG.N,'periodic');
elseif CFG.win == "hamming"
    WIN = hamming(CFG.N,'periodic');
else
    WIN = ones(CFG.N,1);
end
WIN = WIN(:);

% filtro cheby2 bandpass (zero-phase via filtfilt depois)
assert(mod(CFG.ford,2)==0, 'CFG.ford precisa ser par.');
Ncheb = CFG.ford/2;
Rs = 30; % stopband ripple (dB)
[bBP,aBP] = cheby2(Ncheb, Rs, [CFG.wi CFG.wf]/(FS_TARGET/2), 'bandpass');

% resposta do filtro (plot "filtro")
[Hf, Wf] = freqz(bBP, aBP, 4096, FS_TARGET);

figure('Name','Filtro BP (cheby2) - magnitude/fase','Color','w');
subplot(2,1,1);
plot(Wf, 20*log10(abs(Hf)+eps), 'LineWidth',1.2); grid on;
xlim([0 FS_TARGET/2]);
xlabel('Hz'); ylabel('|H| (dB)');
title(sprintf('Cheby2 BP | ord_final=%d (n=%d) | [%.2f %.2f] Hz | Fs=%.1f', CFG.ford, Ncheb, CFG.wi, CFG.wf, FS_TARGET));
subplot(2,1,2);
plot(Wf, unwrap(angle(Hf)), 'LineWidth',1.2); grid on;
xlim([0 FS_TARGET/2]);
xlabel('Hz'); ylabel('fase (rad)');

%% ===================== RUN =====================
Corrs = []; MAEs = []; RMSEs = []; used = 0;

fprintf('\n==================== RUN (BP_THEN_GATE) ====================\n');

for d = 1:numel(DATASETS)
    fprintf('\n[DATASET] %s\n', DATASETS(d).name);

    % --- radar ---
    A = readmatrix(DATASETS(d).radar);
    tpuro = A(:,COL_T_MS);
    xraw  = A(:,COL_HEART);

    tpuro = tpuro(:);
    xraw  = xraw(:);

    t0 = (tpuro - tpuro(1)) / 1000;
    [t0, ia] = unique(t0,'stable');
    xraw = xraw(ia);

    % ======= PRE (sem wrap) =======
    ok = isfinite(t0) & isfinite(xraw);
    t0 = t0(ok); xraw = xraw(ok);

    if numel(xraw) < 8
        fprintf('  [SKIP] poucos pontos raw.\n');
        continue;
    end

    x0 = double(xraw);
    x0 = x0 - mean(x0,'omitnan');

    % --- segmentação por buracos ---
    dt0 = diff(t0);
    brk = find(dt0 > GAP_THR_SEC);
    if isempty(brk)
        segments = [1 numel(t0)];
    else
        segments = [[1; brk+1], [brk; numel(t0)]];
    end

    % --- polar ---
    hasPolar = true;
    try
        fid = fopen(DATASETS(d).polar,'r');
        if fid==-1, error('open fail'); end
        rawL = textscan(fid,'%s','Delimiter','\n');
        fclose(fid);
        rawL = rawL{1};
        if numel(rawL)<2, error('polar vazio'); end

        lines = rawL(2:end);
        ts_str = strings(numel(lines),1);
        HR_polar = nan(numel(lines),1);
        for i=1:numel(lines)
            L = string(lines{i});
            if contains(L,";"), parts = split(L,";"); else, parts = split(L,","); end
            if numel(parts)>=2
                ts_str(i) = strtrim(parts(1));
                HR_polar(i)= str2double(strtrim(parts(2)));
            end
        end

        try
            t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS");
        catch
            t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd HH:mm:ss.SSS");
        end
        t_polar = seconds(t_dt - t_dt(1));
        t_polar = t_polar(:); HR_polar = HR_polar(:);

        mp = isfinite(t_polar) & isfinite(HR_polar);
        t_polar = t_polar(mp); HR_polar = HR_polar(mp);
        [t_polar, iu] = unique(t_polar,'stable');
        HR_polar = HR_polar(iu);

        hasPolar = ~isempty(t_polar);
    catch
        hasPolar = false;
        t_polar = []; HR_polar = [];
    end

    if ~hasPolar
        fprintf('  [SKIP] Polar indisponível.\n');
        continue;
    end

    % ===== RESAMPLE por segmento -> concat =====
    dt = 1/FS_TARGET;
    t_u = [];
    x_u = [];
    for k = 1:size(segments,1)
        i0 = segments(k,1); i1 = segments(k,2);
        ts = t0(i0:i1); xs = x0(i0:i1);

        [ts, iu] = unique(ts,'stable'); xs = xs(iu);
        if numel(ts) < 4, continue; end

        tnew = (ts(1):dt:ts(end))';
        xnew = interp1(ts, xs, tnew, 'linear');

        ok2 = isfinite(tnew) & isfinite(xnew);
        t_u = [t_u; tnew(ok2)];
        x_u = [x_u; xnew(ok2)];
    end
    [t_u, iu] = unique(t_u,'stable');
    x_u = x_u(iu);

    if numel(x_u) < max(MIN_POINTS_AFTER_RESAMPLE, CFG.N)
        fprintf('  [SKIP] Poucos pontos após resample: %d\n', numel(x_u));
        continue;
    end

    % ===== FILTRO (BP) =====
    % evita filtfilt crash por padlen
    padlen = 3*(max(length(aBP), length(bBP)) - 1);
    if numel(x_u) <= padlen
        x_bp = x_u;
    else
        x_bp = filtfilt(bBP,aBP,x_u);
    end
    x_bp = x_bp(:);

    % ===== TRACK (BP_THEN_GATE) + RECON x_g =====
    debug_gate_this = DEBUG_GATE_PLOT && (DATASETS(d).name == DEBUG_GATE_DATASET);

    nfft = CFG.N;
    half = floor(nfft/2) + 1;

    k_min = ceil(F_MIN_HZ * nfft / FS_TARGET) + 1;
    k_min = max(2, k_min);

    hr_min_hz = HR_MIN_BPM/60;
    hr_max_hz = HR_MAX_BPM/60;
    k_hr1 = max(k_min, ceil(hr_min_hz * nfft / FS_TARGET) + 1);
    k_hr2 = min(half,  floor(hr_max_hz * nfft / FS_TARGET) + 1);
    if k_hr2 <= k_hr1
        fprintf('  [SKIP] faixa HR inválida em bins.\n');
        continue;
    end

    gate_pow_ratio = 10^(DELTA_DB/10);

    if DEBUG_RECON_ON
        x_g  = zeros(size(x_bp));
        wsum = zeros(size(x_bp));
    else
        x_g = [];
    end

    t_frames = [];
    f_hz     = [];

    idx0 = 1;
    L = numel(x_bp);
    frame_count = 0;

    while (idx0 + CFG.N - 1) <= L
        frame_count = frame_count + 1;

        xw = x_bp(idx0:idx0+CFG.N-1) .* WIN;
        xw = xw - mean(xw,'omitnan');

        X  = fft(xw, nfft);
        X1 = X(1:half);

        % --- potência pra estimar pico in-band ---
        P1 = abs(X1).^2;
        if SMOOTH_BINS > 1
            P1s = movmean(P1, SMOOTH_BINS, 'Endpoints','shrink');
        else
            P1s = P1;
        end

        % --- gate atua em X1 (magnitude) ---
        X1g = X1;

        if GATE_ON
            Ppk_in = max(P1s(k_hr1:k_hr2));
            if isfinite(Ppk_in) && Ppk_in > 0
                capPow = Ppk_in * gate_pow_ratio;
                if isfinite(GATE_FLOOR) && GATE_FLOOR > 0
                    capPow = max(capPow, GATE_FLOOR);
                end
                capMag = sqrt(capPow);

                Mag1 = abs(X1);
                outband = true(half,1);
                outband(k_hr1:k_hr2) = false;
                outband(1) = false;

                mask_over = outband & (Mag1 > capMag);

                if GATE_MODE == "ZERO"
                    X1g(mask_over) = 0;
                else
                    X1g(mask_over) = X1(mask_over) .* (capMag ./ Mag1(mask_over));
                end

                if debug_gate_this && (frame_count == DEBUG_GATE_FRAME)
                    f = (0:half-1)' * (FS_TARGET/nfft);
                    figure('Name','DEBUG GATE (1 frame)','Color','w');
                    subplot(2,1,1);
                    plot(f, abs(X1), 'LineWidth',1.1); hold on;
                    plot(f, abs(X1g),'LineWidth',1.1);
                    yline(capMag,'--','capMag');
                    xline((k_hr1-1)*FS_TARGET/nfft,'--','HRmin');
                    xline((k_hr2-1)*FS_TARGET/nfft,'--','HRmax');
                    grid on; xlabel('Hz'); ylabel('|X|');
                    legend('|X1|','|X1g|','Location','best');
                    title(sprintf('%s frame %d | capMag=%.4g | kHR=[%d..%d]', DATASETS(d).name, frame_count, capMag, k_hr1, k_hr2));

                    subplot(2,1,2);
                    stem(f, double(mask_over), '.'); grid on;
                    xlabel('Hz'); ylabel('mask\_over');
                    title('Bins fora da faixa HR que estouraram o teto');
                end
            end
        end

        % --- pico no espectro gated (single-sided) ---
        Pg = abs(X1g).^2;

        [~, rel] = max(Pg(k_min:end));
        k0 = k_min + rel - 1;

        % --- harmonic guard (sub-harmônicos) ---
        if CFG.harm_on
            P0 = Pg(k0);
            if isfinite(P0) && P0 > 0
                for m = 2:CFG.harm_max
                    k_est = round(k0 / m);
                    if k_est < k_min, continue; end
                    lo = max(k_min, k_est - CFG.harm_tol);
                    hi = min(numel(Pg), k_est + CFG.harm_tol);
                    [Psub, idx] = max(Pg(lo:hi));
                    ksub = lo + idx - 1;
                    if isfinite(Psub) && (Psub >= CFG.harm_ratio * P0)
                        k0 = ksub;
                        P0 = Psub;
                    end
                end
            end
        end

        % --- interp parabólica em log(P) ---
        if k0 <= 1 || k0 >= numel(Pg)
            delta = 0;
        else
            la = log(Pg(k0-1) + 1e-20);
            lb = log(Pg(k0)   + 1e-20);
            lc = log(Pg(k0+1) + 1e-20);
            denom = (la - 2*lb + lc);
            if denom==0 || ~isfinite(denom)
                delta = 0;
            else
                delta = 0.5*(la - lc)/denom;
                if ~isfinite(delta), delta = 0; end
                delta = max(-0.75, min(0.75, delta));
            end
        end

        f_est = (k0 - 1 + delta) * FS_TARGET / nfft;

        t_center = t_u(idx0) + 0.5*(CFG.N-1)/FS_TARGET;
        t_frames(end+1,1) = t_center; %#ok<AGROW>
        f_hz(end+1,1)     = f_est;    %#ok<AGROW>

        % --- recon overlap-add (debug) ---
        if DEBUG_RECON_ON
            Xfull = zeros(nfft,1);
            Xfull(1:half) = X1g;
            Xfull(half+1:end) = conj(X1g(half-1:-1:2));

            yg = real(ifft(Xfull, nfft));

            x_g(idx0:idx0+CFG.N-1) = x_g(idx0:idx0+CFG.N-1) + yg .* WIN;
            wsum(idx0:idx0+CFG.N-1)= wsum(idx0:idx0+CFG.N-1)+ (WIN.^2);
        end

        idx0 = idx0 + CFG.hop;
    end

    if isempty(t_frames)
        fprintf('  [SKIP] Track vazio.\n');
        continue;
    end

    if DEBUG_RECON_ON
        m = wsum > 0;
        x_g(m) = x_g(m) ./ wsum(m);
    end

    hr_bpm = 60 * f_hz;

    % ===== PLOT SINAIS (pre, bp, bp+gate) =====
    figure('Name', "Sinais - "+DATASETS(d).name, 'Color','w');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    nexttile; plot(t_u, x_u, 'LineWidth',1.0); grid on;
    xlabel('t (s)'); ylabel('amp'); title(DATASETS(d).name + " | PRE (resample)");

    nexttile; plot(t_u, x_bp,'LineWidth',1.0); grid on;
    xlabel('t (s)'); ylabel('amp'); title(DATASETS(d).name + " | BP (filtfilt cheby2)");

    nexttile;
    if DEBUG_RECON_ON
        plot(t_u, x_g,'LineWidth',1.0);
        title(DATASETS(d).name + " | BP + GATE (recon overlap-add)");
    else
        plot(t_u, x_bp,'LineWidth',1.0);
        title(DATASETS(d).name + " | BP + GATE (sem recon) (mostrando BP)");
    end
    grid on; xlabel('t (s)'); ylabel('amp');

    % ===== GRID COMUM + SUAVIZAÇÃO + MÉTRICAS =====
    t_start  = max(min(t_frames), min(t_polar));
    t_end    = min(max(t_frames), max(t_polar));
    if ~(isfinite(t_start) && isfinite(t_end) && t_end > t_start)
        fprintf('  [SKIP] Janela comum inválida.\n');
        continue;
    end

    t_common = (t_start : DT_GRID : t_end)';

    % interp radar -> grid + smooth
    HR_radar_i = interp1(t_frames, hr_bpm, t_common, 'linear', nan);
    win_n = max(1, round(SMOOTH_HR_SEC/DT_GRID));
    HR_radar_grid = movmean(HR_radar_i, win_n, 'omitnan');

    % polar -> grid
    HR_polar_grid = interp1(t_polar, HR_polar, t_common, 'linear', nan);

    m = isfinite(HR_radar_grid) & isfinite(HR_polar_grid);
    n = sum(m);
    if n < MIN_METRIC_SAMPLES
        fprintf('  [SKIP] Métrica inválida (n=%d).\n', n);
        continue;
    end

    e = HR_radar_grid(m) - HR_polar_grid(m);
    MAE  = mean(abs(e));
    RMSE = sqrt(mean(e.^2));
    R    = corr(HR_radar_grid(m), HR_polar_grid(m));

    used = used + 1;
    Corrs(end+1,1) = R;    %#ok<AGROW>
    MAEs(end+1,1)  = MAE;  %#ok<AGROW>
    RMSEs(end+1,1) = RMSE; %#ok<AGROW>

    fprintf('  OK | corr=%.4f | MAE=%.3f | RMSE=%.3f | n=%d\n', R, MAE, RMSE, n);

    figure('Name', "HR (grid) - "+DATASETS(d).name, 'Color','w');
    plot(t_common, HR_polar_grid, 'LineWidth',1.2); hold on;
    plot(t_common, HR_radar_grid, 'LineWidth',1.2);
    grid on;
    xlabel('t (s)'); ylabel('HR (bpm)');
    legend('Polar (grid)','Radar (grid)','Location','best');
    title(sprintf('%s | corr=%.3f MAE=%.2f RMSE=%.2f', DATASETS(d).name, R, MAE, RMSE));

    % ===== deixar no workspace (último dataset que rodou OK) =====
    assignin('base','t_u',t_u);
    assignin('base','x_u',x_u);
    assignin('base','x_f',x_bp);
    if DEBUG_RECON_ON, assignin('base','x_g',x_g); end
    assignin('base','t_frames',t_frames);
    assignin('base','hr_bpm',hr_bpm);
end

fprintf('\n==================== SUMMARY ====================\n');
if used == 0
    fprintf('Nenhum dataset válido.\n');
else
    fprintf('Datasets usados: %d/%d\n', used, numel(DATASETS));
    fprintf('corr_mean=%.4f | MAE_mean=%.3f | RMSE_mean=%.3f\n', mean(Corrs), mean(MAEs), mean(RMSEs));
end
