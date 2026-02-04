clc; clear; close all;

%% ===================== PARÂMETROS =====================
INPUT_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\datasets\Figshare\dta\GDN0008\GDN0008_5_TiltDown.mat';
INPUT_MODE = 'radar';   % 'radar' | 'ecg1' | 'ecg2'

wi = 0.5;      % Hz (CWT/PSD display)
wf = 3.0;      % Hz

SRATE_TARGET = 25;  % Hz

% --- GATE espectral ---
GATE_ON     = 1;
HR_MIN_BPM  = 55;
HR_MAX_BPM  = 210;
DELTA_DB    = -3;   % thr = picoHR + DELTA_DB  (negativo = mais agressivo)
SMOOTH_BINS = 3;    % suaviza Pdb e W
GATE_FLOOR  = 0.0;  % 0 = mata, >0 atenua

% --- Bandpass (após gate) ---
BP_ORDER = 6;
BP_WI    = 0.7;     % Hz
BP_WF    = 2.5;     % Hz

% --- CWT ---
VOICES_PER_OCT = 48;
WIN_SCALO_SEC  = 50;     % abas de 50s
DB_FLOOR_PLOT  = -25;    % dB rel pico global

% --- Ridge (área) ---
BAND_W_HZ     = 0.15;  % Hz
DB_FLOOR_BAND = -25;   % dB rel pico local (piso p/ área)

% --- ECG gabarito (só pra avaliação/plot) ---
ECG_MIN_DIST_SEC = 0.35;
ECG_MIN_PROM     = 2.0;
HR_MIN_BPM_ECG   = 30;
HR_MAX_BPM_ECG   = 220;

% --- Punish / evento ---
RIDGE_LAMBDA        = 0.25;
F_JUMP_HZ           = 0.15;
BURNIN_SEC          = 0;
K_EVENT             = 30;
PUNISH_CHAIN_GAMMA  = 0.8;
PUNISH_CHAIN_MAX    = 12;

% --- Comparador harmônico (f/2) ---
HARM_ON        = 1;      % 1 = habilita
HARM_REL_DB    = -8;     % aceita f/2 se área(f/2) >= área(f)*10^(HARM_REL_DB/10)
HARM_MIN_BPM   = 40;     % sanity range (evita jogar pra respiração sem querer)
HARM_MAX_BPM   = 220;

%% ===================== LEITURA =====================
[tpuro_ms, heart, baseName] = load_phase_for_algorithm(INPUT_PATH);

heart = heart(:);
heart = heart - mean(heart,'omitnan');

%% ===================== SRATE ROBUSTO =====================
t0 = (tpuro_ms - tpuro_ms(1))/1000;
t0 = t0(:);
if numel(t0) < 2, error('Sinal muito curto.'); end

dt_pos = diff(t0); dt_pos = dt_pos(dt_pos > 0);
if isempty(dt_pos), error('Tempo inválido: diff(t) sem valores positivos.'); end

dt0    = median(dt_pos);
srate0 = 1/dt0;

%% ===================== RESAMPLE =====================
t0dur = (numel(heart)-1)/srate0;
[p,q] = rat(SRATE_TARGET/srate0, 1e-12);

heart = resample(heart, p, q);
srate = srate0 * p / q;
t     = (0:numel(heart)-1)'/srate;
if t(end) > t0dur, t(end) = t0dur; end

% validações
if wf >= (srate/2), error('wf precisa ser < Nyquist. wf=%.3f | Nyq=%.3f', wf, srate/2); end
if wi <= 0 || wi >= wf, error('wi/wf inválidos: wi=%.3f, wf=%.3f', wi, wf); end

heart_puro = heart(:);

%% ===================== GATE (ANTES DO BP) =====================
if GATE_ON
    [heart_gate, gate_dbg] = gate_by_power_hr( ...
        heart_puro, srate, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);
else
    heart_gate = heart_puro;
    gate_dbg = struct('bpm',[],'Pdb_before',[],'Pdb_after',[],'Pdb_hr_peak',NaN,'Pdb_thr',NaN);
end

%% ===================== BANDPASS (SOS estável) =====================
if BP_WF >= (srate/2)
    error('BP_WF precisa ser < Nyquist. BP_WF=%.3f | Nyquist=%.3f', BP_WF, srate/2);
end
if BP_WI <= 0 || BP_WI >= BP_WF
    error('BP_WI/BP_WF inválidos: BP_WI=%.3f, BP_WF=%.3f', BP_WI, BP_WF);
end

Wn = [BP_WI BP_WF]/(srate/2);

use_sos_direct = false;
try
    [sosBP, gBP] = butter(BP_ORDER, Wn, 'bandpass', 'sos');
    use_sos_direct = true;
catch
    use_sos_direct = false;
end

if ~use_sos_direct
    [z,pz,kz] = butter(BP_ORDER, Wn, 'bandpass');
    sosBP = zp2sos(z,pz,kz);
    gBP   = 1;
end

if exist('sosfiltfilt','file') == 2
    heart = sosfiltfilt(sosBP, gBP, heart_gate);
else
    [bBP,aBP] = sos2tf(sosBP, gBP);
    heart = filtfilt(bBP, aBP, heart_gate);
end

%% ===================== PSD overlay =====================
[fpsd, P_puro] = psd_one_sided(heart_puro, srate);
[~,    P_gate] = psd_one_sided(heart_gate, srate);
[~,    P_bp]   = psd_one_sided(heart,      srate);

mpsd = (fpsd >= wi) & (fpsd <= wf);

figure('Name','PSD overlay (a.u.)','Color','w');
plot(fpsd(mpsd), P_puro(mpsd), 'LineWidth', 1.2); hold on;
plot(fpsd(mpsd), P_gate(mpsd), 'LineWidth', 1.2);
plot(fpsd(mpsd), P_bp(mpsd),   'LineWidth', 1.6);
grid on; xlabel('Frequência (Hz)'); ylabel('Potência (a.u.)');
title(sprintf('PSD | %s | Fs=%.2f | wi=%.2f wf=%.2f | GATE=%d', baseName, srate, wi, wf, GATE_ON));
legend('puro (pós-resample)','pós-gate','pós-gate+BP','Location','best');

figure('Name','PSD overlay (dB)','Color','w');
plot(fpsd(mpsd), 10*log10(P_puro(mpsd)+eps), 'LineWidth', 1.2); hold on;
plot(fpsd(mpsd), 10*log10(P_gate(mpsd)+eps), 'LineWidth', 1.2);
plot(fpsd(mpsd), 10*log10(P_bp(mpsd)+eps),   'LineWidth', 1.6);
grid on; xlabel('Frequência (Hz)'); ylabel('Potência (dB a.u.)');
title(sprintf('PSD (dB) | %s | Fs=%.2f | GATE=%d', baseName, srate, GATE_ON));
legend('puro (pós-resample)','pós-gate','pós-gate+BP','Location','best');

%% ===================== CWT =====================
[wt, f_bin] = cwt(heart,'amor',srate, ...
    'FrequencyLimits',[wi wf], ...
    'VoicesPerOctave',VOICES_PER_OCT);

f_bin = f_bin(:);
if numel(f_bin) > 1 && f_bin(2) < f_bin(1)
    f_bin = flipud(f_bin);
    wt    = flipud(wt);
end

% potência (usar single ajuda memória/velocidade)
P  = single(abs(wt).^2);     % nb x nt
clear wt;

nb = numel(f_bin);
nt = size(P,2);

df = gradient(f_bin); df(df<=0) = eps;

%% ===================== PRECOMPUTES RIDGE =====================
burnin_frames = max(1, round(BURNIN_SEC * srate));

% end_idx(i): último bin j com f(j) <= f(i)+BAND_W_HZ
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

% mapa "metade da freq" para um bin inicial viável (aprox nearest)
half_idx_float = interp1(f_bin, 1:nb, f_bin/2, 'linear', NaN);
half_start_idx = int32(round(half_idx_float));
half_start_idx(~isfinite(half_idx_float)) = int32(0);
half_start_idx(half_start_idx < 1 | half_start_idx > nb) = int32(0);

FLOOR_FRAC = 10^(DB_FLOOR_BAND/10);
HARM_THR_LIN = 10^(HARM_REL_DB/10);

%% ===================== RIDGE (ÁREA + PUNISH + HARM f/2) =====================
imax    = zeros(nt,1,'int32');
band_lo = zeros(nt,1,'int32');
band_hi = zeros(nt,1,'int32');

prev_f_normal  = NaN;
prev_was_event = false;

sorted_peak = NaN(nt,1,'single');
n_sorted    = 0;

punish_chain = 0;

% buffers para evitar re-alocação dentro do loop
kpeak = zeros(nb,1,'int32');
fpeak = zeros(nb,1);

for tt = 1:nt
    pcol = double(P(:,tt));      % double aqui pra estabilidade dos cálculos
    pmax = max(pcol) + eps;

    if n_sorted == 0
        base_p = pmax;
    else
        base_p = double(median_sorted(sorted_peak, n_sorted));
    end

    is_event = (pmax > K_EVENT * base_p);

    % área acima do piso relativo (pico local)
    pnorm      = pcol / pmax;
    excess_lin = max(pnorm - FLOOR_FRAC, 0);

    c = [0; cumsum(excess_lin .* df)];  % c tem nb+1
    ei = double(end_idx);              % pra indexar no c
    score_area = c(ei+1) - c(1:nb);
    scoreN = score_area / (max(score_area) + eps);

    % pico local por faixa (versão leve: loop curto, janela pequena)
    for i = 1:nb
        e = double(end_idx(i));
        [~, rel] = max(pcol(i:e));
        kpeak(i) = int32(i + rel - 1);
    end
    fpeak(:) = f_bin(double(kpeak));

    [~, i_best_area] = max(scoreN);

    score_total = scoreN;

    lambda_eff = RIDGE_LAMBDA;
    if is_event, lambda_eff = 0; end

    prev_for_punish = prev_f_normal;
    if (~is_event) && prev_was_event
        prev_for_punish = NaN;
    end

    punish_applied = false;
    if (tt > burnin_frames) && ~isnan(prev_for_punish) && (lambda_eff > 0)
        k = min(PUNISH_CHAIN_MAX, punish_chain);
        lambda_chain = lambda_eff * exp(-PUNISH_CHAIN_GAMMA * k);

        jump = (fpeak - prev_for_punish) / (F_JUMP_HZ + eps);
        score_total = score_total - lambda_chain * (jump.^2);
        punish_applied = true;
    end

    [~, i_best] = max(score_total);

    % --- HARMÔNICO: testa se f/2 é viável (corrige quando caiu no 2º harmônico) ---
    if HARM_ON
        % pega a freq escolhida (pelo pico do melhor candidato)
        k_main = kpeak(i_best);
        i_half = half_start_idx(k_main);

        if i_half > 0
            % critério por área (e/ou score_total), com sanity de BPM
            area_main = score_area(i_best) + eps;
            area_half = score_area(double(i_half)) + eps;

            f_half_peak = f_bin(double(kpeak(double(i_half))));
            bpm_half = 60*f_half_peak;

            if isfinite(bpm_half) && (bpm_half >= HARM_MIN_BPM) && (bpm_half <= HARM_MAX_BPM)
                if (area_half >= HARM_THR_LIN * area_main)
                    % troca para f/2
                    i_best = double(i_half);
                end
            end
        end
    end

    band_lo(tt) = int32(i_best);
    band_hi(tt) = end_idx(i_best);
    imax(tt)    = kpeak(i_best);

    p_peak = pcol(double(imax(tt)));

    if ~is_event
        prev_f_normal = f_bin(double(imax(tt)));
        [sorted_peak, n_sorted] = sorted_insert(sorted_peak, n_sorted, single(p_peak));
    end

    if punish_applied && (~is_event) && (i_best ~= i_best_area)
        punish_chain = min(PUNISH_CHAIN_MAX, punish_chain + 1);
    else
        punish_chain = 0;
    end

    prev_was_event = is_event;
end

freq_hz  = f_bin(double(imax));
freq_bpm = 60*freq_hz;

%% ===================== SCALOGRAMA (dB) =====================
Pmax = double(max(P(:))) + eps;
Sdb  = 10*log10(double(P)/Pmax + eps);
Sdb(Sdb < DB_FLOOR_PLOT) = NaN;

%% ===================== GABARITO ECG + MÉTRICAS =====================
[hr_ref_bpm, ~] = hr_gabarito_from_ecg(INPUT_PATH, t, ECG_MIN_DIST_SEC, ECG_MIN_PROM, HR_MIN_BPM_ECG, HR_MAX_BPM_ECG);

if isempty(hr_ref_bpm)
    warning('Sem gabarito ECG neste arquivo. Vai plotar só o ridge.');
else
    ok = isfinite(freq_bpm) & isfinite(hr_ref_bpm);
    est = freq_bpm(ok);
    ref = hr_ref_bpm(ok);

    err  = est - ref;
    RMSE = sqrt(mean(err.^2));
    MAE  = mean(abs(err));
    if numel(est) >= 3
        C = corrcoef(est, ref);
        CORR = C(1,2);
    else
        CORR = NaN;
    end

    fprintf('\n=== MÉTRICAS (ridge vs ECG) ===\n');
    fprintf('RMSE=%.3f | MAE=%.3f | CORR=%.4f\n', RMSE, MAE, CORR);
end

%% ===================== MOVMEAN 7s + PLOT =====================
WIN_MEAN_SEC = 7;
winS = max(1, round(WIN_MEAN_SEC * srate));

freq_bpm_smooth = movmean(freq_bpm, winS, 'Endpoints','shrink');

if ~isempty(hr_ref_bpm)
    hr_ref_bpm_smooth = movmean(hr_ref_bpm, winS, 'Endpoints','shrink');
else
    hr_ref_bpm_smooth = [];
end

if ~isempty(hr_ref_bpm)
    okR = isfinite(freq_bpm) & isfinite(hr_ref_bpm);
    okS = isfinite(freq_bpm_smooth) & isfinite(hr_ref_bpm_smooth);

    estR = freq_bpm(okR);          refR = hr_ref_bpm(okR);
    estS = freq_bpm_smooth(okS);   refS = hr_ref_bpm_smooth(okS);

    RMSE_raw = sqrt(mean((estR-refR).^2));
    MAE_raw  = mean(abs(estR-refR));
    CORR_raw = corr(estR, refR, 'Rows','complete');

    RMSE_sm  = sqrt(mean((estS-refS).^2));
    MAE_sm   = mean(abs(estS-refS));
    CORR_sm  = corr(estS, refS, 'Rows','complete');

    fprintf('\n=== MOVMEAN %.1fs (win=%d) ===\n', WIN_MEAN_SEC, winS);
    fprintf('RAW   -> RMSE=%.3f | MAE=%.3f | CORR=%.4f\n', RMSE_raw, MAE_raw, CORR_raw);
    fprintf('SMOOTH-> RMSE=%.3f | MAE=%.3f | CORR=%.4f\n', RMSE_sm,  MAE_sm,  CORR_sm);
end

% >>> sem plot do ECG "raw" (somente ECG smooth, se existir)
figure('Name', sprintf('HR: movmean %.1fs vs ECG', WIN_MEAN_SEC), 'Color','w');
plot(t, freq_bpm_smooth, 'LineWidth', 2.0); hold on;

if ~isempty(hr_ref_bpm_smooth)
    plot(t, hr_ref_bpm_smooth, '--', 'LineWidth', 2.0);
    legend(sprintf('Radar movmean %.1fs', WIN_MEAN_SEC), ...
           sprintf('ECG movmean %.1fs', WIN_MEAN_SEC), ...
           'Location','best');
else
    legend(sprintf('Radar movmean %.1fs', WIN_MEAN_SEC), 'Location','best');
end

grid on; xlabel('t (s)'); ylabel('HR (bpm)');
title(sprintf('%s | movmean %.1fs | GATE=%d', baseName, WIN_MEAN_SEC, GATE_ON));

%% ===================== PLOT EM ABAS: scalograma + ridge (+ gabarito smooth) =====================
win_N = max(1, round(WIN_SCALO_SEC * srate));
nT    = size(P,2);
nWin  = ceil(nT / win_N);

fig = figure('Name', sprintf('%s | CWT + Ridge | abas %ds', baseName, WIN_SCALO_SEC), ...
             'NumberTitle','off', 'Color','w');
tg = uitabgroup(fig);

for k = 1:nWin
    i1 = (k-1)*win_N + 1;
    i2 = min(k*win_N, nT);

    tab = uitab(tg, 'Title', sprintf('%02d: %.0f-%.0fs', k, t(i1), t(i2)));

    ax1 = axes('Parent', tab); ax1.Color = 'w';
    Sseg  = Sdb(:, i1:i2);
    maskA = isfinite(Sseg);

    hImg = imagesc(ax1, t(i1:i2), f_bin, Sseg);
    set(hImg, 'AlphaData', maskA);
    axis(ax1,'xy'); grid(ax1,'on');
    xlabel(ax1,'t (s)'); ylabel(ax1,'Freq (Hz)');
    title(ax1, sprintf('Scalograma (dB) + Ridge | Janela %d (%.0f–%.0f s)', k, t(i1), t(i2)));
    cb = colorbar(ax1); ylabel(cb,'dB (rel.)');
    caxis(ax1, [DB_FLOOR_PLOT 0]);

    drawnow;
    pos1 = ax1.Position;

    ax2 = axes('Parent', tab, ...
        'Position', pos1, ...
        'Color','none', ...
        'YAxisLocation','right', ...
        'XAxisLocation','bottom', ...
        'Box','off');
    hold(ax2,'on');

    plot(ax2, t(i1:i2), freq_bpm(i1:i2), 'k', 'LineWidth', 2.6, 'HandleVisibility','off');
    plot(ax2, t(i1:i2), freq_bpm(i1:i2), 'w', 'LineWidth', 1.3, 'DisplayName','Radar ridge');

    if ~isempty(hr_ref_bpm_smooth)
        plot(ax2, t(i1:i2), hr_ref_bpm_smooth(i1:i2), 'c--', 'LineWidth', 2.0, ...
            'DisplayName',sprintf('ECG movmean %.1fs', WIN_MEAN_SEC));
    end

    ax2.XLim = [t(i1) t(i2)];
    ax2.YLim = [f_bin(1) f_bin(end)]*60;
    ax2.XTick = ax1.XTick;
    ylabel(ax2,'HR (bpm)');
    ax2.YColor = [1 1 1];

    legend(ax2, 'Location','northeast');
end

%% ===================== FUNÇÕES =====================
function [tpuro_ms, heart, baseName] = load_phase_for_algorithm(path_in)
    [~, baseName, ext] = fileparts(path_in);
    ext = lower(ext);

    if ~strcmp(ext,'.mat') && ~strcmp(ext,'.csv')
        error('Entrada inválida: use .csv ou .mat.');
    end

    if strcmp(ext,'.csv')
        data = readmatrix(path_in);
        tpuro_ms = data(:,1);
        phase    = data(:,3);
        heart    = wrap_phase(phase);
        return;
    end

    S = load(path_in);

    if evalin('base','exist(''INPUT_MODE'',''var'')')
        mode = evalin('base','INPUT_MODE');
    else
        mode = 'radar';
    end
    mode = lower(string(mode));

    switch mode
        case "ecg1"
            if ~(isfield(S,'tfm_ecg1') && isfield(S,'fs_ecg'))
                error('INPUT_MODE=ecg1, mas faltam tfm_ecg1/fs_ecg no MAT.');
            end
            x  = double(S.tfm_ecg1(:));
            fs = double(S.fs_ecg);

        case "ecg2"
            if ~(isfield(S,'tfm_ecg2') && isfield(S,'fs_ecg'))
                error('INPUT_MODE=ecg2, mas faltam tfm_ecg2/fs_ecg no MAT.');
            end
            x  = double(S.tfm_ecg2(:));
            fs = double(S.fs_ecg);

        otherwise % radar
            if ~(isfield(S,'radar_i') && isfield(S,'radar_q') && isfield(S,'fs_radar'))
                error('INPUT_MODE=radar, mas faltam radar_i/radar_q/fs_radar no MAT.');
            end
            I  = double(S.radar_i(:));
            Q  = double(S.radar_q(:));
            fs = double(S.fs_radar);

            phase = unwrap(atan2(Q, I));
            heart = wrap_phase(phase);

            tpuro_ms = ((0:numel(I)-1)'/fs) * 1000;
            return;
    end

    tpuro_ms = ((0:numel(x)-1)'/fs) * 1000;
    heart = wrap_phase(unwrap(x));
end

function ph_wrapped = wrap_phase(ph)
    ph = ph(:);
    ok = isfinite(ph);
    if ~all(ok)
        ph(~ok) = interp1(find(ok), ph(ok), find(~ok), 'linear', 'extrap');
    end
    ph_wrapped = mod(ph + pi, 2*pi) - pi;
    ph_wrapped = ph_wrapped - mean(ph_wrapped,'omitnan');
end

function [f, P1] = psd_one_sided(x, Fs)
    x = x(:);
    x = x - mean(x,'omitnan');
    N = numel(x);
    NFFT = 2^nextpow2(N);

    X = fft(x, NFFT);
    P2 = (abs(X)/N).^2;

    P1 = P2(1:NFFT/2+1);
    if numel(P1) > 2
        P1(2:end-1) = 2*P1(2:end-1);
    end

    f = (0:NFFT/2)' * (Fs/NFFT);
end

%% ===== GATE: broad (Pdb suavizado) + spike (Pdb bruto) =====
function [xatt, dbg] = gate_by_power_hr(x, FS, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR)
    x = x(:);
    N = numel(x);

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
    if ~any(mHR)
        xatt = x;
        dbg = struct('bpm',[],'Pdb_before',[],'Pdb_after',[],'Pdb_hr_peak',NaN,'Pdb_thr',NaN);
        return;
    end

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

    kk  = 0:floor(N/2);
    ff  = kk*(FS/N);
    bpmU = 60*ff;

    P_before = abs(X(1:numel(kk))).^2;
    P_after  = abs(Xg(1:numel(kk))).^2;

    if N > 1 && numel(P_before) > 2
        P_before(2:end-1) = 2*P_before(2:end-1);
        P_after(2:end-1)  = 2*P_after(2:end-1);
    end

    dbg.bpm = bpmU(:);
    dbg.Pdb_before = 10*log10(P_before/Pref + eps);
    dbg.Pdb_after  = 10*log10(P_after /Pref + eps);

    mHRu = (dbg.bpm >= HR_MIN_BPM) & (dbg.bpm <= HR_MAX_BPM);
    if any(mHRu)
        dbg.Pdb_hr_peak = max(dbg.Pdb_before(mHRu), [], 'omitnan');
    else
        dbg.Pdb_hr_peak = NaN;
    end
    dbg.Pdb_thr = dbg.Pdb_hr_peak + DELTA_DB;
end

function [hr_ref_bpm, dbg] = hr_gabarito_from_ecg(path_in, t_radar, minDistSec, minProm, hrMin, hrMax)
    dbg = struct('n_peaks',0);

    if ~endsWith(lower(path_in), '.mat')
        hr_ref_bpm = [];
        return;
    end

    S = load(path_in);
    if ~(isfield(S,'tfm_ecg1') && isfield(S,'fs_ecg'))
        hr_ref_bpm = [];
        return;
    end

    ecg = double(S.tfm_ecg1(:));
    fs  = double(S.fs_ecg);

    if isempty(ecg) || ~isfinite(fs) || fs <= 0
        hr_ref_bpm = [];
        return;
    end

    w_baseline = max(3, round(0.20 * fs));
    ecg0 = ecg - movmedian(ecg, w_baseline, 'Endpoints','shrink');

    f1 = 5;  f2 = 25;
    if f2 >= fs/2, f2 = 0.9*(fs/2); end
    if f1 <= 0 || f1 >= f2
        hr_ref_bpm = NaN(size(t_radar));
        return;
    end

    [b,a] = butter(2, [f1 f2]/(fs/2), 'bandpass');
    ecg_f = filtfilt(b,a, ecg0);
    ecg_f = ecg_f / (mad(ecg_f,1) + eps);

    d = [0; diff(ecg_f)];
    d = d / (mad(d,1) + eps);
    mwi = movmean(d.^2, max(3, round(0.12*fs)), 'Endpoints','shrink');

    thr = median(mwi) + 3.0*mad(mwi,1);
    minDistSec2 = max(minDistSec, 60/hrMax);
    minDistSmp  = max(1, round(minDistSec2 * fs));
    prom_adapt  = max(minProm, 1.5*mad(mwi,1));

    [~, locs_env] = findpeaks(mwi, ...
        'MinPeakDistance',  minDistSmp, ...
        'MinPeakHeight',    thr, ...
        'MinPeakProminence',prom_adapt);

    dbg.n_peaks = numel(locs_env);
    if dbg.n_peaks < 3
        hr_ref_bpm = NaN(size(t_radar));
        return;
    end

    w_ref = max(1, round(0.08 * fs));
    locs_R = zeros(size(locs_env));
    n = numel(ecg_f);

    for k = 1:numel(locs_env)
        c = locs_env(k);
        i1 = max(1, c - w_ref);
        i2 = min(n, c + w_ref);
        [~, im] = max(abs(ecg_f(i1:i2)));
        locs_R(k) = i1 + im - 1;
    end

    locs_R = unique(locs_R, 'stable');
    dbg.n_peaks = numel(locs_R);
    if dbg.n_peaks < 3
        hr_ref_bpm = NaN(size(t_radar));
        return;
    end

    t_ecg = (0:numel(ecg_f)-1)'/fs;
    t_R   = t_ecg(locs_R);

    RR = diff(t_R);
    if isempty(RR)
        hr_ref_bpm = NaN(size(t_radar));
        return;
    end

    hr_inst = 60 ./ RR;
    t_hr = t_R(1:end-1) + RR/2;

    ok = isfinite(hr_inst) & isfinite(t_hr) & (hr_inst >= hrMin) & (hr_inst <= hrMax);

    RR_ok = RR(ok);
    if numel(RR_ok) >= 5
        medRR = median(RR_ok);
        ok = ok & (RR > 0.5*medRR) & (RR < 1.8*medRR);
    end

    t_hr    = t_hr(ok);
    hr_inst = hr_inst(ok);

    if numel(hr_inst) < 3
        hr_ref_bpm = NaN(size(t_radar));
        return;
    end

    hr_ref_bpm = interp1(t_hr, hr_inst, t_radar, 'linear', 'extrap');
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
