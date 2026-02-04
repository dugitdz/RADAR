clc; clear; close all;

%% ===================== PATHS =====================
CSV_PATH      = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\phase.csv';
POLAR_HR_PATH = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\POLARH10.txt';
TES_PATH      = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes.csv';   % formato: t_sec,HR

%% ===================== CSV COLUNAS =====================
COL_T_MS  = 1;
COL_PHASE = 2;

%% ===================== CONFIG =====================
gap_thr     = 0.5;     % s
MIN_SEG_SEC = 2.0;     % s

FS_TARGET = 25.0;      % Hz
IMETH     = 'linear';

WRAP_ON   = 0;         % 1: wrap+mean | 0: mean apenas

FLOW = 'BP_THEN_GATE'; % 'BP_THEN_GATE' ou 'GATE_THEN_BP'

% GATE
GATE_ON      = 1;
HR_MIN_BPM   = 65;
HR_MAX_BPM   = 210;
DELTA_DB     = -5;
SMOOTH_BINS  = 3;
GATE_FLOOR   = 0.0;

% BP (BUTTER SOS)
BP_ON     = 1;
BP_ORDER  = 2;
BP_WI     = 0.3;   % Hz  (SEM CLAMP. usa isto mesmo.)
BP_WF     = 2.5;       % Hz

% CWT
WI       = 0.6;        % Hz
WF       = 3.0;        % Hz
VOICES   = 48;
CWT_WAVE = 'amor';

% RIDGE
BAND_W_HZ      = 0.15;
DB_FLOOR_BAND  = -25;
RIDGE_LAMBDA   = 0.3;
F_JUMP_HZ      = 0.15;
BURNIN_SEC     = 5;
K_EVENT        = 30;
PCHAIN_GAMMA   = 0.4;
PCHAIN_MAX     = 6;

% HARM
HARM_ON      = 1;
HARM_REL_DB  = -6;
HARM_MIN_BPM = 65;
HARM_MAX_BPM = 220;

% suavização final
MOVMEAN_SEC = 7.0;

%% ===================== LEITURA CSV =====================
A = readmatrix(CSV_PATH);
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

%% ===================== SEGMENTAÇÃO POR BURACOS =====================
segments = segment_by_gaps(t0, gap_thr);

%% ===================== LOOP: RESAMPLE -> (FLOW) -> CWT -> RIDGE =====================
dt = 1/FS_TARGET;

tseg_list  = {};
hrseg_list = {};

best_len = -inf;
spec = struct('t',[],'pure',[],'gate',[],'flow',[]);

for s = 1:size(segments,1)
    i0 = segments(s,1); i1 = segments(s,2);

    ts = t0(i0:i1);
    xs = x0(i0:i1);

    [ts, iu] = unique(ts, 'stable');
    xs = xs(iu);

    if numel(ts) < 8, continue; end

    tnew = (ts(1):dt:ts(end))';
    xnew = interp1(ts, xs, tnew, IMETH);

    ok = isfinite(tnew) & isfinite(xnew);
    tnew = tnew(ok);
    xnew = xnew(ok);

    if numel(tnew) < 16, continue; end
    if (tnew(end) - tnew(1)) < MIN_SEG_SEC, continue; end

    xnew = force_finite_vector(xnew);

    % =================== SINAIS PARA PLOT (COERENTE COM FLOW) ===================
    x_pure_spec = xnew;

    if GATE_ON == 1
        x_gate_spec = gate_by_power_hr_full(x_pure_spec, FS_TARGET, ...
            HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
    else
        x_gate_spec = x_pure_spec;
    end

    x_flow_spec = x_pure_spec;

    if strcmpi(FLOW,'BP_THEN_GATE')
        if BP_ON == 1
            x_flow_spec = apply_bp_sos_unitgain(x_flow_spec, FS_TARGET, BP_ORDER, BP_WI, BP_WF);
        end
        if GATE_ON == 1
            x_flow_spec = gate_by_power_hr_full(x_flow_spec, FS_TARGET, ...
                HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
        end
    else
        % ======== GATE_THEN_BP ========
        if GATE_ON == 1
            x_flow_spec = gate_by_power_hr_full(x_flow_spec, FS_TARGET, ...
                HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
        end
        if BP_ON == 1
            % >>> SEM CLAMP. usa BP_WI/BP_WF direto.
            x_flow_spec = apply_bp_sos_unitgain(x_flow_spec, FS_TARGET, BP_ORDER, BP_WI, BP_WF);
        end
    end

    seg_len = tnew(end) - tnew(1);
    if seg_len > best_len
        best_len = seg_len;
        spec.t    = tnew(:);
        spec.pure = x_pure_spec(:);
        spec.gate = x_gate_spec(:);
        spec.flow = x_flow_spec(:);
    end

    % =================== PIPELINE FINAL (HR) ===================
    x_fin = xnew;

    if (GATE_ON ~= 1) && (BP_ON ~= 1)
        % nada
    elseif GATE_ON ~= 1
        if BP_ON == 1
            x_fin = apply_bp_sos_unitgain(x_fin, FS_TARGET, BP_ORDER, BP_WI, BP_WF);
        end
    elseif BP_ON ~= 1
        x_fin = gate_by_power_hr_full(x_fin, FS_TARGET, ...
            HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
    else
        if strcmpi(FLOW,'BP_THEN_GATE')
            x_fin = apply_bp_sos_unitgain(x_fin, FS_TARGET, BP_ORDER, BP_WI, BP_WF);
            x_fin = gate_by_power_hr_full(x_fin, FS_TARGET, ...
                HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
        else
            x_fin = gate_by_power_hr_full(x_fin, FS_TARGET, ...
                HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
            % >>> SEM CLAMP
            x_fin = apply_bp_sos_unitgain(x_fin, FS_TARGET, BP_ORDER, BP_WI, BP_WF);
        end
    end

    x_fin = force_finite_vector(x_fin);

    % ---------- CWT ----------
    [wt, fbin] = cwt(x_fin, CWT_WAVE, FS_TARGET, ...
        'FrequencyLimits',[WI WF], ...
        'VoicesPerOctave', VOICES);

    fbin = fbin(:);
    if numel(fbin)>1 && fbin(2)<fbin(1)
        fbin = flipud(fbin);
        wt   = flipud(wt);
    end

    P = single(abs(wt).^2);
    clear wt;
    P(~isfinite(P)) = 0;

    % ---------- RIDGE ----------
    freq_hz = ridge_area_punish_event_harm(P, fbin, FS_TARGET, ...
        BAND_W_HZ, DB_FLOOR_BAND, RIDGE_LAMBDA, F_JUMP_HZ, BURNIN_SEC, ...
        K_EVENT, PCHAIN_GAMMA, PCHAIN_MAX, ...
        HARM_ON, HARM_REL_DB, HARM_MIN_BPM, HARM_MAX_BPM);

    hrseg = 60*freq_hz(:);

    if MOVMEAN_SEC > 0
        Wseg = max(1, round(MOVMEAN_SEC * FS_TARGET));
        hrseg = movmean(hrseg, Wseg, 'Endpoints','shrink');
    end

    tseg_list{end+1}  = tnew(:);
    hrseg_list{end+1} = hrseg(:);
end

%% ===================== POLAR =====================
[t_polar, HR_polar] = read_txt_polar_flex(POLAR_HR_PATH);
tG = t_polar(:);
hP = HR_polar(:);

%% ===================== TES =====================
[tT, hT] = read_tes_csv(TES_PATH);

%% ===================== LIGAR PONTAS (RADAR) =====================
[t_all, h_all] = concat_segments(tseg_list, hrseg_list);
hC_onP = interp1(t_all, h_all, tG, 'linear', NaN);
hC_onT = interp1(t_all, h_all, tT, 'linear', NaN);
hP_onT = interp1(tG, hP, tT, 'linear', NaN);

%% ===================== MÉTRICAS: RADAR vs POLAR =====================
mE = isfinite(hC_onP) & isfinite(hP);
e  = hC_onP(mE) - hP(mE);

RMSE = sqrt(mean(e.^2,'omitnan'));
MAE  = mean(abs(e),'omitnan');

if nnz(mE) >= 5
    CORR = corr(hC_onP(mE), hP(mE), 'Rows','complete');
else
    CORR = NaN;
end

fprintf('\n===== Polar vs CWT@Polar =====\n');
fprintf('N Polar = %d | N válido = %d\n', numel(tG), nnz(mE));
fprintf('RMSE=%.4f | MAE=%.4f | Corr=%.4f\n', RMSE, MAE, CORR);

%% ===================== MÉTRICAS EXTRA =====================
mRT = isfinite(hC_onT) & isfinite(hT);
eRT = hC_onT(mRT) - hT(mRT);
RMSE_RT = sqrt(mean(eRT.^2,'omitnan'));
MAE_RT  = mean(abs(eRT),'omitnan');
if nnz(mRT) >= 5
    CORR_RT = corr(hC_onT(mRT), hT(mRT), 'Rows','complete');
else
    CORR_RT = NaN;
end

mTP = isfinite(hT) & isfinite(hP_onT);
eTP = hT(mTP) - hP_onT(mTP);
RMSE_TP = sqrt(mean(eTP.^2,'omitnan'));
MAE_TP  = mean(abs(eTP),'omitnan');
if nnz(mTP) >= 5
    CORR_TP = corr(hT(mTP), hP_onT(mTP), 'Rows','complete');
else
    CORR_TP = NaN;
end

fprintf('\n===== EXTRA =====\n');
fprintf('RADAR vs TES:   N=%d | RMSE=%.4f | MAE=%.4f | Corr=%.4f\n', nnz(mRT), RMSE_RT, MAE_RT, CORR_RT);
fprintf('TES vs POLAR:   N=%d | RMSE=%.4f | MAE=%.4f | Corr=%.4f\n', nnz(mTP), RMSE_TP, MAE_TP, CORR_TP);

%% ===================== PLOTS EM ABAS =====================
fig = figure('Color','w');
tg = uitabgroup(fig);

tab1 = uitab(tg, 'Title','Curvas (tempo)');
ax1 = axes(tab1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, tG, hP, 'c--', 'LineWidth', 1.6);
plot(ax1, tG, hC_onP, 'k-', 'LineWidth', 1.4);
plot(ax1, tT, hT, 'm-', 'LineWidth', 1.2);
xlabel(ax1,'t (s)'); ylabel(ax1,'HR (bpm)');
title(ax1, sprintf(['POLAR vs RADAR vs TES | ' ...
    'RADAR-POLAR: MAE=%.3f RMSE=%.3f Corr=%.3f | ' ...
    'RADAR-TES: MAE=%.3f RMSE=%.3f Corr=%.3f | ' ...
    'TES-POLAR: MAE=%.3f RMSE=%.3f Corr=%.3f'], ...
    MAE, RMSE, CORR, MAE_RT, RMSE_RT, CORR_RT, MAE_TP, RMSE_TP, CORR_TP));
legend(ax1,'POLAR','RADAR (interp@POLAR)','TES','Location','best');

tab2 = uitab(tg, 'Title','Potência (bins) radar');
ax2 = axes(tab2); hold(ax2,'on'); grid(ax2,'on');

% >>>>>> pedido: tirar RAW do plot e forçar 0 dB no pico mais alto
if ~isempty(spec.t) && numel(spec.gate) >= 16 && numel(spec.flow) >= 16
    [bpm, Pdb_gate, Pdb_flow] = power_bins_overlay2( ...
        spec.gate, spec.flow, FS_TARGET, 1);

    plot(ax2, bpm, Pdb_gate, 'LineWidth', 1.2);
    plot(ax2, bpm, Pdb_flow, 'LineWidth', 1.2);

    % DEBUG opcional: mostra onde o FLOW "subiu" vs gate (não depende do 0 dB do topo)
    % plot(ax2, bpm, (Pdb_flow - Pdb_gate), 'LineWidth', 1.2);

    xlim(ax2, [0 300]);
    xlabel(ax2,'BPM (bins FFT)');
    ylabel(ax2,'Potência (dB rel. ao maior pico)');
    title(ax2, sprintf('Maior segmento: %.1fs | Fs=%.1f | FLOW=%s (0 dB no pico)', best_len, FS_TARGET, FLOW));
    legend(ax2,'Gate','Saída FLOW','Location','best');
else
    text(ax2, 0.05, 0.5, 'Sem segmento válido para FFT/potência.', 'Units','normalized');
    axis(ax2,'off');
end

%% ===================== FUNÇÕES =====================

function segments = segment_by_gaps(t, gap_thr)
    t = t(:);
    n = numel(t);
    if n < 2, segments = [1 n]; return; end
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
        if numel(tS) < 2, continue; end
        [tS, iu] = unique(tS, 'stable');
        hS = hS(iu);
        t_all = [t_all; tS]; %#ok<AGROW>
        h_all = [h_all; hS]; %#ok<AGROW>
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
    if ~any(isfinite(x)), x(:) = 0; end
end

% ======= BP com ganho unitário (evita "subir" bins) =======
function x = apply_bp_sos_unitgain(x, Fs, ord, wi, wf)
    x = force_finite_vector(x);

    if ~isfinite(Fs) || Fs<=0, return; end
    if wf >= Fs/2, return; end
    if wi <= 0 || wi >= wf, return; end

    Wn = [wi wf]/(Fs/2);

    try
        [sosBP, gBP] = butter(ord, Wn, 'bandpass', 'sos');
    catch
        [z,p,k] = butter(ord, Wn, 'bandpass');
        sosBP = zp2sos(z,p,k);
        gBP = 1;
    end

    % Normaliza ganho em f0 (centro geométrico)
    f0 = sqrt(wi*wf);
    w0 = 2*pi*(f0/Fs);

    try
        H0 = freqz(sosBP, gBP, w0);
        m0 = abs(H0);
        if isfinite(m0) && m0 > 0
            gBP = gBP / m0;
        end
    catch
        % segue sem normalização se freqz falhar
    end

    if exist('sosfiltfilt','file') == 2
        x = sosfiltfilt(sosBP, gBP, x);
    else
        [b,a] = sos2tf(sosBP, gBP);
        x = filtfilt(b,a,x);
    end

    x = force_finite_vector(x);
end

function xatt = gate_by_power_hr_full(x, FS, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR)
    x = force_finite_vector(x);
    N = numel(x);
    if N < 8 || ~isfinite(FS) || FS<=0
        xatt = x; return;
    end

    SPIKE_HALF  = 1;
    DO_W_SMOOTH = (SMOOTH_BINS > 1);

    X = fft(x);

    k = (0:N-1).';
    f = k*(FS/N);
    f_abs = min(f, FS - f);
    bpm_full = 60*f_abs;

    P2   = abs(X).^2;
    Pref = max(P2) + eps;

    Pdb_raw = 10*log10(P2/Pref + eps);
    if SMOOTH_BINS > 1
        Pdb_s = movmean(Pdb_raw, SMOOTH_BINS, 'omitnan');
    else
        Pdb_s = Pdb_raw;
    end

    mHR = (bpm_full >= HR_MIN_BPM) & (bpm_full <= HR_MAX_BPM);
    if ~any(mHR), xatt = x; return; end

    Pdb_hr_peak = max(Pdb_s(mHR), [], 'omitnan');
    thr = Pdb_hr_peak + DELTA_DB;

    mOUT = ~mHR;
    mask_broad = mOUT & (Pdb_s   > thr);
    mask_spike = mOUT & (Pdb_raw > thr);

    if SPIKE_HALF > 0 && any(mask_spike)
        ker = ones(2*SPIKE_HALF+1,1);
        mask_spike = conv(double(mask_spike), ker, 'same') > 0;
        mask_spike = logical(mask_spike);
    end

    mask = mask_broad | mask_spike;

    W = ones(N,1);
    W(mask) = GATE_FLOOR;
    W(mHR)  = 1;

    if DO_W_SMOOTH
        W = movmean(W, SMOOTH_BINS, 'omitnan');
        W = min(1, max(GATE_FLOOR, W));
        W(mHR) = 1;
        W(mask_spike) = GATE_FLOOR;
    end

    Xg = X .* W;
    xatt = real(ifft(Xg, 'symmetric'));
    xatt = force_finite_vector(xatt);
end

% ======= POTÊNCIA (2 curvas) | 0 dB no pico do FLOW =======
function [bpm, Pdb_gate, Pdb_flow] = power_bins_overlay2(x_gate, x_flow, Fs, smooth_bins)
    x_gate = force_finite_vector(x_gate);
    x_flow = force_finite_vector(x_flow);

    N = min([numel(x_gate), numel(x_flow)]);
    x_gate = x_gate(1:N);
    x_flow = x_flow(1:N);

    Xg = fft(x_gate);
    Xf = fft(x_flow);

    Pg = abs(Xg).^2;
    Pf = abs(Xf).^2;

    k = (0:N-1).';
    f = k*(Fs/N);
    f_abs = min(f, Fs - f);
    bpm_full = 60*f_abs;

    nh  = floor(N/2) + 1;
    bpm = bpm_full(1:nh);

    Pg = Pg(1:nh);
    Pf = Pf(1:nh);

    % referência = pico do FLOW (antes de suavizar)
    Pref = max(Pf) + eps;

    Pdb_gate = 10*log10(Pg/Pref + eps);
    Pdb_flow = 10*log10(Pf/Pref + eps);

    if smooth_bins > 1
        Pdb_gate = movmean(Pdb_gate, smooth_bins, 'omitnan');
        Pdb_flow = movmean(Pdb_flow, smooth_bins, 'omitnan');
    end

    % garante FLOW com topo exatamente em 0 dB (mesmo após smooth)
    peak_flow = max(Pdb_flow);
    if isfinite(peak_flow)
        Pdb_gate = Pdb_gate - peak_flow;
        Pdb_flow = Pdb_flow - peak_flow;
    end
end

% ===== ridge (igual ao seu) =====
function freq_hz = ridge_area_punish_event_harm(P, f_bin, srate, BAND_W_HZ, DB_FLOOR_BAND, RIDGE_LAMBDA, F_JUMP_HZ, BURNIN_SEC, K_EVENT, GAMMA, KMAX, ...
                                                HARM_ON, HARM_REL_DB, HARM_MIN_BPM, HARM_MAX_BPM)
    f_bin = double(f_bin(:));
    nb = numel(f_bin);
    nt = size(P,2);

    df = gradient(f_bin);
    df(df<=0) = eps;

    burnin_frames = max(1, round(BURNIN_SEC * srate));

    end_idx = zeros(nb,1,'int32');
    j = 1;
    for i = 1:nb
        if j < i, j = i; end
        f_end = f_bin(i) + BAND_W_HZ;
        while (j < nb) && (f_bin(j+1) <= f_end)
            j = j + 1;
        end
        end_idx(i) = int32(j);
    end

    half_start_idx = int32(zeros(nb,1));
    if HARM_ON == 1
        half_idx_float = interp1(f_bin, 1:nb, f_bin/2, 'linear', NaN);
        tmp = int32(round(half_idx_float));
        tmp(~isfinite(half_idx_float)) = int32(0);
        tmp(tmp < 1 | tmp > nb) = int32(0);
        half_start_idx = tmp;
    end

    FLOOR_FRAC   = 10^(DB_FLOOR_BAND/10);
    HARM_THR_LIN = 10^(HARM_REL_DB/10);

    imax = zeros(nt,1,'int32');
    prev_f_normal  = NaN;
    prev_was_event = false;

    sorted_peak = NaN(nt,1,'single');
    n_sorted    = 0;

    punish_chain = 0;

    kpeak = zeros(nb,1,'int32');
    fpeak = zeros(nb,1);

    score_area = zeros(nb,1);
    scoreN     = zeros(nb,1);
    score_total= zeros(nb,1);

    for tt = 1:nt
        pcol = double(P(:,tt));
        pmax = max(pcol) + eps;

        if n_sorted == 0
            base_p = pmax;
        else
            base_p = double(median_sorted(sorted_peak, n_sorted));
        end

        is_event = (pmax > K_EVENT * base_p);

        pnorm      = pcol / pmax;
        excess_lin = max(pnorm - FLOOR_FRAC, 0);

        c = [0; cumsum(excess_lin .* df)];
        ei = double(end_idx);
        score_area(:) = c(ei+1) - c(1:nb);
        scoreN(:)     = score_area / (max(score_area) + eps);

        for i = 1:nb
            e = double(end_idx(i));
            [~, rel] = max(pcol(i:e));
            kpeak(i) = int32(i + rel - 1);
        end
        fpeak(:) = f_bin(double(kpeak));

        [~, i_best_area] = max(scoreN);
        score_total(:) = scoreN;

        lambda_eff = RIDGE_LAMBDA;
        if is_event, lambda_eff = 0; end

        prev_for_punish = prev_f_normal;
        if (~is_event) && prev_was_event
            prev_for_punish = NaN;
        end

        punish_applied = false;
        if (tt > burnin_frames) && ~isnan(prev_for_punish) && (lambda_eff > 0)
            k = min(KMAX, punish_chain);
            lambda_chain = lambda_eff * exp(-GAMMA * k);
            jump = (fpeak - prev_for_punish) / (F_JUMP_HZ + eps);
            score_total(:) = score_total(:) - lambda_chain * (jump.^2);
            punish_applied = true;
        end

        [~, i_best] = max(score_total);

        if HARM_ON == 1
            k_main = kpeak(i_best);
            i_half = half_start_idx(k_main);

            if i_half > 0
                area_main = score_area(i_best) + eps;
                area_half = score_area(double(i_half)) + eps;

                f_half_peak = f_bin(double(kpeak(double(i_half))));
                bpm_half = 60*f_half_peak;

                if isfinite(bpm_half) && (bpm_half >= HARM_MIN_BPM) && (bpm_half <= HARM_MAX_BPM)
                    if (area_half >= HARM_THR_LIN * area_main)
                        i_best = double(i_half);
                    end
                end
            end
        end

        imax(tt) = kpeak(i_best);

        p_peak = pcol(double(imax(tt)));
        if ~is_event
            prev_f_normal = f_bin(double(imax(tt)));
            [sorted_peak, n_sorted] = sorted_insert(sorted_peak, n_sorted, single(p_peak));
        end

        if punish_applied && (~is_event) && (i_best ~= i_best_area)
            punish_chain = min(KMAX, punish_chain + 1);
        else
            punish_chain = 0;
        end

        prev_was_event = is_event;
    end

    freq_hz = f_bin(double(imax));
end

function m = median_sorted(a, n)
    if n <= 0, m = NaN; return; end
    if mod(n,2)==1
        m = a((n+1)/2);
    else
        m = 0.5*(a(n/2) + a(n/2 + 1));
    end
end

function [a, n] = sorted_insert(a, n, x)
    if ~isfinite(x), return; end
    if n == 0
        a(1) = x; n = 1; return;
    end
    lo = 1; hi = n;
    while lo <= hi
        mid = floor((lo+hi)/2);
        if x < a(mid), hi = mid - 1; else, lo = mid + 1; end
    end
    idx = lo;
    if idx <= n
        a(idx+1:n+1) = a(idx:n);
    end
    a(idx) = x;
    n = n + 1;
end

function [t_sec, HR] = read_txt_polar_flex(p)
    fid = fopen(p,'r');
    if fid == -1
        t_sec = []; HR = []; return;
    end
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
        t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS");
    catch
        t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd HH:mm:ss.SSS");
    end

    t_sec = seconds(t_dt - t_dt(1));
    t_sec = t_sec(:); HR = HR(:);

    m = isfinite(t_sec) & isfinite(HR);
    t_sec = t_sec(m); HR = HR(m);

    [t_sec, iu] = unique(t_sec,'stable');
    HR = HR(iu);
end

function [tT, hT] = read_tes_csv(path)
    M = readmatrix(path);
    if isempty(M) || size(M,2) < 2
        tT = []; hT = [];
        return;
    end
    tT = double(M(:,1));
    hT = double(M(:,2));

    tT = tT(:); hT = hT(:);
    ok = isfinite(tT) & isfinite(hT);
    tT = tT(ok); hT = hT(ok);

    [tT, iu] = unique(tT, 'stable');
    hT = hT(iu);
end
