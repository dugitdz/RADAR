clc; clear; close all;

%% ===================== PATHS =====================
CSV_PATH      = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv';
POLAR_HR_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH3.txt';
TES_PATH      = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes.csv';   % formato: t_sec,HR

%% ===================== CSV COLUNAS =====================
COL_T_MS  = 1;   % tempo em ms
COL_PHASE = 4;   % fase (rad) já pronta

%% ============================================================
%  AUTO-APLICAR: melhor config por CORRELAÇÃO (do teu GRID)
%  - Usa:
%       grid_out/grid_aggregate_by_config.csv
%       grid_out/grid_all_results.mat  (contém CFG)
%% ============================================================
USE_BEST_CORR_FROM_GRID = true;

GRID_OUT_DIR  = fullfile(pwd, 'grid_out');
AG_CSV_PATH   = fullfile(GRID_OUT_DIR, 'grid_aggregate_by_config.csv');
GRID_MAT_PATH = fullfile(GRID_OUT_DIR, 'grid_all_results.mat');

% Se quiser forçar um id específico, set para número; senão, NaN:
FORCE_CFG_ID = 643;  % NaN para não forçar

best_cfg_id = NaN;

if USE_BEST_CORR_FROM_GRID
    if ~exist(AG_CSV_PATH,'file'),   error('Não achei: %s', AG_CSV_PATH); end
    if ~exist(GRID_MAT_PATH,'file'), error('Não achei: %s', GRID_MAT_PATH); end

    AG = readtable(AG_CSV_PATH);

    ok = isfinite(AG.CORRg) & isfinite(AG.ONE_MINUS_CORRg) & (AG.Nds_ok >= 1);
    AG2 = AG(ok,:);
    if isempty(AG2), error('AG não tem configs válidas com CORRg finita.'); end

    [~, ii] = max(AG2.CORRg);
    best_cfg_id = AG2.cfg_id(ii);

    if isfinite(FORCE_CFG_ID)
        best_cfg_id = FORCE_CFG_ID;
    end

    S = load(GRID_MAT_PATH, 'CFG');
    CFG = S.CFG;

    if best_cfg_id < 1 || best_cfg_id > numel(CFG)
        error('best_cfg_id=%d fora do range do CFG (1..%d).', best_cfg_id, numel(CFG));
    end

    p = CFG(best_cfg_id);

    fprintf('\n[GRID] cfg_id=%d | CORRg_max=%.5f | ONE_MINUS=%.5f | Nds_ok=%d\n', ...
        best_cfg_id, AG2.CORRg(ii), AG2.ONE_MINUS_CORRg(ii), AG2.Nds_ok(ii));

    % ---------------- PARAMETROS DO GRID ----------------
    gap_thr     = p.gap_thr;
    MIN_SEG_SEC = p.MIN_SEG_SEC;

    BP_ON     = p.BP_ON;
    FILT_TYPE = string(p.FILT_TYPE);
    ORD_FINAL = p.ORD_FINAL;
    WI        = p.WI;
    WF        = p.WF;

    GATE_ON     = p.GATE_ON;
    HR_MIN_BPM  = p.HR_MIN_BPM;
    HR_MAX_BPM  = p.HR_MAX_BPM;
    DELTA_DB    = p.DELTA_DB;
    SMOOTH_BINS = p.SMOOTH_BINS;
    GATE_FLOOR  = p.GATE_FLOOR;

    wi             = p.wi;
    wf             = p.wf;
    VOICES_PER_OCT = 16;

    MOVMEAN_SEC = 7;
    WRAP_ON     = p.WRAP_ON;

else
    % ---------------- MANUAL ----------------
    gap_thr     = 0.5;
    MIN_SEG_SEC = 2.0;

    BP_ON      = 1;
    FILT_TYPE  = "butter";
    ORD_FINAL  = 4;
    WI         = 0.5;
    WF         = 2.5;

    GATE_ON      = 1;
    HR_MIN_BPM   = 70;
    HR_MAX_BPM   = 210;
    DELTA_DB     = 0;
    SMOOTH_BINS  = 9;
    GATE_FLOOR   = 0.0;

    wi = 0.0; wf = 3.5; VOICES_PER_OCT = 48;

    MOVMEAN_SEC = 7.0;
    WRAP_ON     = 1;
end

%% ===================== Resample / Plot =====================
FS_TARGET     = 25.0;      % Hz
IMETH         = 'linear';
DB_FLOOR_PLOT = -40;       % dB

%% ===================== LEITURA CSV (RADAR) =====================
A = readmatrix(CSV_PATH);
tpuro_ms = A(:, COL_T_MS);
phase0   = A(:, COL_PHASE);

tpuro_ms = tpuro_ms(:);
phase0   = phase0(:);

t0 = (tpuro_ms - tpuro_ms(1))/1000; % s
[t0, ia] = unique(t0(:), 'stable');
phase0   = phase0(ia);

% pré-tratamento
if WRAP_ON == 1
    x0 = wrap_phase(phase0);
else
    x0 = double(phase0(:));
    okx = isfinite(x0);
    if ~all(okx)
        x0(~okx) = interp1(find(okx), x0(okx), find(~okx), 'linear', 'extrap');
    end
    x0 = x0 - mean(x0,'omitnan');
end

%% ===================== SEGMENTAÇÃO POR BURACOS =====================
segments = segment_by_gaps(t0, gap_thr);

%% ===================== FILTRO FINAL (COEFS) =====================
if BP_ON == 1
    [bBP, aBP, minLenFilt] = make_final_bp(FILT_TYPE, ORD_FINAL, WI, WF, FS_TARGET);
else
    bBP = []; aBP = []; minLenFilt = 0;
end

%% ===================== LOOP: RESAMPLE -> GATE -> (BP) -> CWT =====================
dt = 1/FS_TARGET;

t_plot = [];
P_all  = [];
f_bin  = [];

% overlay no scalograma (com NaN no gap)
t_hr_vis  = [];
hr_vis_sm = [];

% curva por segmento (pra ligar depois)
tseg_list     = {};
hrseg_list_sm = {};

% debug (1o segmento válido)
dbg_saved = false;
dbg_t = []; dbg_x_before = []; dbg_x_after = [];
dbg_bpm = []; dbg_Pdb_before = []; dbg_Pdb_after = [];
dbg_hr_peak_db = NaN; dbg_thr_db = NaN;

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

    x_pre = xnew;

    % (1) GATE
    if GATE_ON == 1
        [x_gate, gate_dbg] = gate_by_power_hr( ...
            x_pre, FS_TARGET, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR);

        if ~dbg_saved
            dbg_saved    = true;
            dbg_t        = tnew(:);
            dbg_x_before = x_pre(:);
            dbg_x_after  = x_gate(:);

            dbg_bpm        = gate_dbg.bpm(:);
            dbg_Pdb_before = gate_dbg.Pdb_before(:);
            dbg_Pdb_after  = gate_dbg.Pdb_after(:);
            dbg_hr_peak_db = gate_dbg.Pdb_hr_peak;
            dbg_thr_db     = gate_dbg.Pdb_thr;
        end
    else
        x_gate = x_pre;
    end

    % (2) BP final
    if BP_ON == 1
        if numel(x_gate) < minLenFilt, continue; end
        x_for_cwt = filtfilt(bBP, aBP, x_gate);
    else
        x_for_cwt = x_gate;
    end

    % (3) CWT
    [wt, f] = cwt(x_for_cwt, 'morse', FS_TARGET, ...
        'VoicesPerOctave', VOICES_PER_OCT, ...
        'FrequencyLimits', [wi wf]);

    f = f(:);
    if numel(f) > 1 && f(2) < f(1)
        f  = flipud(f);
        wt = flipud(wt);
    end

    Pseg = abs(wt).^2;
    Pseg(~isfinite(Pseg)) = 0;

    % HR por coluna (imax)
    [~, imax] = max(Pseg, [], 1);
    hrseg = 60 * f(imax);
    tseg  = tnew(:);
    hrseg = hrseg(:);

    % movmean só dentro do segmento
    hrseg_sm = hrseg;
    if MOVMEAN_SEC > 0
        Wseg = max(1, round(MOVMEAN_SEC * FS_TARGET));
        if numel(hrseg_sm) >= Wseg
            hrseg_sm = movmean(hrseg_sm, Wseg, 'omitnan');
        end
    end

    % guarda
    tseg_list{end+1}     = tseg;
    hrseg_list_sm{end+1} = hrseg_sm;

    % ---- scalograma global ----
    if isempty(f_bin)
        f_bin  = f;
        P_all  = Pseg;
        t_plot = tnew;

        t_hr_vis  = tseg;
        hr_vis_sm = hrseg_sm;
    else
        if numel(f) ~= numel(f_bin) || any(abs(f - f_bin) > 1e-9)
            Pseg2 = interp1(f, Pseg, f_bin, 'linear', 'extrap');
            Pseg2(~isfinite(Pseg2)) = 0;
            Pseg2(Pseg2 < 0) = 0;
        else
            Pseg2 = Pseg;
        end

        P_all  = [P_all, nan(size(P_all,1), 1), Pseg2]; %#ok<AGROW>
        t_plot = [t_plot; (t_plot(end)+dt); tnew];     %#ok<AGROW>

        t_hr_vis  = [t_hr_vis; (t_hr_vis(end)+dt); tseg]; %#ok<AGROW>
        hr_vis_sm = [hr_vis_sm; NaN; hrseg_sm];          %#ok<AGROW>
    end
end

%% ===================== CURVA CWT CONTÍNUA (segmentos ligados) =====================
[t_all, h_all] = concat_segments(tseg_list, hrseg_list_sm);

%% ===================== POLAR =====================
[t_polar, HR_polar] = read_txt_polar_flex(POLAR_HR_PATH);
tG = t_polar(:);
hP = HR_polar(:);

%% ===================== TES =====================
[t_tes, HR_tes] = read_csv_time_hr(TES_PATH);
tT = t_tes(:);
hT = HR_tes(:);

%% ===================== AMOSTRAGENS NO TEMPO DAS REFERÊNCIAS =====================
% (mesma curva CWT; só muda a grade de tempo)
hC_onP = interp1(t_all, h_all, tG, 'linear', NaN);
hC_onT = interp1(t_all, h_all, tT, 'linear', NaN);

% TES @ Polar (comparação direta TES vs Polar)
hT_onP = interp1(tT, hT, tG, 'linear', NaN);

%% ===================== MÉTRICAS (CWT vs Polar) =====================
mE = isfinite(hC_onP) & isfinite(hP);
e  = hC_onP(mE) - hP(mE);

RMSE = sqrt(mean(e.^2,'omitnan'));
MAE  = mean(abs(e),'omitnan');
BIAS = mean(e,'omitnan');
SDE  = std(e,'omitnan');
MAXABS = max(abs(e),[],'omitnan');

if nnz(mE) >= 5
    C = corrcoef(hC_onP(mE), hP(mE));
    CORR = C(1,2);
else
    CORR = NaN;
end

fprintf('\n===== ERROS (CWT@Polar - Polar) [segmentos ligados] =====\n');
fprintf('N Polar total = %d | N válido = %d\n', numel(tG), nnz(mE));
fprintf('RMSE=%.4f | MAE=%.4f | BIAS=%.4f | STD=%.4f | max|e|=%.4f | Corr=%.4f\n', ...
    RMSE, MAE, BIAS, SDE, MAXABS, CORR);

%% ===================== MÉTRICAS (CWT vs TES) =====================
mE2 = isfinite(hC_onT) & isfinite(hT);
e2  = hC_onT(mE2) - hT(mE2);

RMSE2 = sqrt(mean(e2.^2,'omitnan'));
MAE2  = mean(abs(e2),'omitnan');
BIAS2 = mean(e2,'omitnan');
SDE2  = std(e2,'omitnan');
MAXABS2 = max(abs(e2),[],'omitnan');

if nnz(mE2) >= 5
    C2 = corrcoef(hC_onT(mE2), hT(mE2));
    CORR2 = C2(1,2);
else
    CORR2 = NaN;
end

fprintf('\n===== ERROS (CWT@TES - TES) [segmentos ligados] =====\n');
fprintf('N TES total = %d | N válido = %d\n', numel(tT), nnz(mE2));
fprintf('RMSE=%.4f | MAE=%.4f | BIAS=%.4f | STD=%.4f | max|e|=%.4f | Corr=%.4f\n', ...
    RMSE2, MAE2, BIAS2, SDE2, MAXABS2, CORR2);

%% ===================== MÉTRICAS (TES@Polar vs Polar) =====================
mTP = isfinite(hT_onP) & isfinite(hP);
eTP = hT_onP(mTP) - hP(mTP);

RMSE_TP = sqrt(mean(eTP.^2,'omitnan'));
MAE_TP  = mean(abs(eTP),'omitnan');
BIAS_TP = mean(eTP,'omitnan');
SDE_TP  = std(eTP,'omitnan');
MAXABS_TP = max(abs(eTP),[],'omitnan');

if nnz(mTP) >= 5
    CTP = corrcoef(hT_onP(mTP), hP(mTP));
    CORR_TP = CTP(1,2);
else
    CORR_TP = NaN;
end

fprintf('\n===== ERROS (TES@Polar - Polar) [TES interpolado em tPolar] =====\n');
fprintf('N Polar total = %d | N válido = %d\n', numel(tG), nnz(mTP));
fprintf('RMSE=%.4f | MAE=%.4f | BIAS=%.4f | STD=%.4f | max|e|=%.4f | Corr=%.4f\n', ...
    RMSE_TP, MAE_TP, BIAS_TP, SDE_TP, MAXABS_TP, CORR_TP);

%% ===================== dB SCALOGRAMA =====================
if isempty(P_all) || isempty(t_plot) || isempty(f_bin)
    Sdb = [];
else
    Pmax = max(P_all(:), [], 'omitnan') + eps;
    Sdb  = 10*log10(P_all/Pmax + eps);
    Sdb(Sdb < DB_FLOOR_PLOT) = NaN;
end

%% ===================== FIG COM ABAS =====================
fig = figure('Name','Radar HR Pipeline (Tabs)','Color','w','Position',[60 60 1550 900]);
tg = uitabgroup(fig);

% ---- Aba 1: CWT + refs ----
tab1 = uitab(tg,'Title','CWT + refs');
ax1 = axes('Parent',tab1); ax1.Color = 'w';

if ~isempty(Sdb)
    hImg = imagesc(ax1, t_plot, f_bin*60, Sdb);
    set(hImg, 'AlphaData', isfinite(Sdb));
    axis(ax1,'xy'); grid(ax1,'on');
    xlabel(ax1,'t (s)'); ylabel(ax1,'Freq (bpm)');

    ttl_gate = "";
    if GATE_ON==1
        ttl_gate = sprintf(' | GATE HR[%d..%d] \\Delta=%gdB', HR_MIN_BPM, HR_MAX_BPM, DELTA_DB);
    end

    if BP_ON == 1
        title(ax1, sprintf('CWT (BP %s ord=%d, %.2f-%.2f Hz) | dB(0..%d) | Fs=%.1f | WRAP=%d%s | MINSEG=%.1fs | cfg=%s', ...
            FILT_TYPE, ORD_FINAL, WI, WF, DB_FLOOR_PLOT, FS_TARGET, WRAP_ON, ttl_gate, MIN_SEG_SEC, string(best_cfg_id)));
    else
        title(ax1, sprintf('CWT (SEM BP) | dB(0..%d) | Fs=%.1f | WRAP=%d%s | MINSEG=%.1fs | cfg=%s', ...
            DB_FLOOR_PLOT, FS_TARGET, WRAP_ON, ttl_gate, MIN_SEG_SEC, string(best_cfg_id)));
    end

    cb = colorbar(ax1); ylabel(cb,'dB (rel.)');
    caxis(ax1,[DB_FLOOR_PLOT 0]);

    hold(ax1,'on');
    hPolarPlot = plot(ax1, tG, hP, 'c--', 'LineWidth', 1.6);
    hTesPlot   = plot(ax1, tT, hT, 'm--', 'LineWidth', 1.6);
    hCwtRidge  = plot(ax1, t_hr_vis, hr_vis_sm, 'w-', 'LineWidth', 1.2);

    legend(ax1, [hPolarPlot hTesPlot hCwtRidge], {'Polar HR','TES HR','HR\_CWT (movmean/seg)'}, 'Location','northeast');
else
    text(ax1,0.1,0.5,'Sem scalograma (P\_all vazio).','Units','normalized');
    axis(ax1,'off');
end

% ---- Aba 2: Comparação geral (1 CWT, sem "dois CWT") ----
tab2 = uitab(tg,'Title','Comparação (geral)');
ax2 = axes('Parent',tab2); ax2.Color = 'w';

plot(ax2, tG, hP, 'c--','LineWidth',1.6); hold(ax2,'on');
plot(ax2, tT, hT, 'm--','LineWidth',1.6);
plot(ax2, t_all, h_all, 'k-','LineWidth',1.4);

grid(ax2,'on'); xlabel(ax2,'t (s)'); ylabel(ax2,'HR (bpm)');
title(ax2, sprintf(['Polar + TES + CWT (segmentos ligados) | movmean(seg)=%.1fs | MINSEG=%.1fs | ' ...
    'CWT@Polar: RMSE=%.3f MAE=%.3f Corr=%.3f | CWT@TES: RMSE=%.3f MAE=%.3f Corr=%.3f | cfg=%s'], ...
    MOVMEAN_SEC, MIN_SEG_SEC, RMSE, MAE, CORR, RMSE2, MAE2, CORR2, string(best_cfg_id)));
legend(ax2,'Polar','TES','CWT contínuo','Location','best');

% ---- Aba 3: TES vs Polar (como você pediu) ----
tab3 = uitab(tg,'Title','TES vs Polar');
ax3TP = axes('Parent',tab3); ax3TP.Color = 'w';

plot(ax3TP, tG, hP, 'c--','LineWidth',1.6); hold(ax3TP,'on');
plot(ax3TP, tG, hT_onP, 'm-','LineWidth',1.2);

grid(ax3TP,'on'); xlabel(ax3TP,'t (s)'); ylabel(ax3TP,'HR (bpm)');
title(ax3TP, sprintf('TES@Polar vs Polar | RMSE=%.3f | MAE=%.3f | Corr=%.3f', RMSE_TP, MAE_TP, CORR_TP));
legend(ax3TP,'Polar','TES interpolado em tPolar','Location','best');

% ---- Aba 4: Debug ----
tab4 = uitab(tg,'Title','Debug');
tlo = tiledlayout(tab4, 2, 1, 'TileSpacing','compact','Padding','compact');

axD1 = nexttile(tlo,1); axD1.Color='w';
if dbg_saved
    plot(axD1, dbg_t, dbg_x_before, 'k-', 'LineWidth', 1.0); hold(axD1,'on');
    plot(axD1, dbg_t, dbg_x_after,  'r-', 'LineWidth', 1.0);
    grid(axD1,'on'); xlabel(axD1,'t (s)'); ylabel(axD1,'amplitude (a.u.)');
    title(axD1, sprintf('Tempo: antes vs atenuado | BP=%d GATE=%d', BP_ON, GATE_ON));
    legend(axD1,'antes','atenuado','Location','best');
else
    text(axD1,0.1,0.5,'Sem segmento válido para debug','Units','normalized'); axis(axD1,'off');
end

axD2 = nexttile(tlo,2); axD2.Color='w';
if dbg_saved && ~isempty(dbg_bpm)
    plot(axD2, dbg_bpm, dbg_Pdb_before, 'k-', 'LineWidth', 1.0); hold(axD2,'on');
    plot(axD2, dbg_bpm, dbg_Pdb_after,  'r-', 'LineWidth', 1.0);

    if isfinite(dbg_hr_peak_db), yline(axD2, dbg_hr_peak_db, 'k--', 'LineWidth', 1.2); end
    if isfinite(dbg_thr_db),     yline(axD2, dbg_thr_db,     'r--', 'LineWidth', 1.2); end

    grid(axD2,'on'); xlabel(axD2,'Frequência (bpm)'); ylabel(axD2,'Potência (dB rel. ao pico do ANTES)');
    title(axD2, sprintf('Espectro dB: antes vs depois | \\Delta=%gdB', DELTA_DB));
    ylim(axD2,[DB_FLOOR_PLOT 0]); xlim(axD2,[0 240]);
    legend(axD2,'antes','atenuado','pico HR','teto','Location','best');
else
    text(axD2,0.1,0.5,'Sem espectro debug','Units','normalized'); axis(axD2,'off');
end

%% ===================== FUNÇÕES (só as usadas) =====================

function x = wrap_phase(ph)
    ph = double(ph(:));
    ok = isfinite(ph);
    if ~all(ok)
        ph(~ok) = interp1(find(ok), ph(ok), find(~ok), 'linear', 'extrap');
    end
    x = mod(ph + pi, 2*pi) - pi;
    x = x - mean(x,'omitnan');
end

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

function [b, a, minLenFilt] = make_final_bp(FILT_TYPE, ORD_FINAL, WI, WF, FS)
    if ORD_FINAL <= 0 || mod(ORD_FINAL,2) ~= 0
        error('ORD_FINAL deve ser par e >0.');
    end
    nyq = FS/2;
    if WI <= 0 || WF >= nyq || WI >= WF
        error('Cortes inválidos: WI=%.3f, WF=%.3f, Nyq=%.3f', WI, WF, nyq);
    end

    switch lower(char(FILT_TYPE))
        case 'butter'
            Wn = [WI WF] / nyq;
            [b, a] = butter(ORD_FINAL/2, Wn, 'bandpass');
        otherwise
            error('FILT_TYPE não suportado: %s', FILT_TYPE);
    end

    nf = max(length(a), length(b)) - 1;
    minLenFilt = 3*nf + 1;
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

function [t_sec, HR] = read_csv_time_hr(p)
    if ~exist(p,'file')
        warning('Não achei TES: %s', p);
        t_sec = []; HR = []; return;
    end

    % tenta: vírgula, depois ;, depois fallback
    try
        M = readmatrix(p, 'Delimiter', ',');
        if isempty(M) || size(M,2) < 2 || all(isnan(M(:)))
            M = readmatrix(p, 'Delimiter', ';');
        end
        if isempty(M) || size(M,2) < 2 || all(isnan(M(:)))
            M = readmatrix(p);
        end
    catch
        M = readmatrix(p);
    end

    if isempty(M) || size(M,2) < 2
        t_sec = []; HR = []; return;
    end

    t_sec = M(:,1);
    HR    = M(:,2);

    t_sec = t_sec(:); HR = HR(:);
    m = isfinite(t_sec) & isfinite(HR);
    t_sec = t_sec(m); HR = HR(m);

    [t_sec, iu] = unique(t_sec, 'stable');
    HR = HR(iu);
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

function [xatt, dbg] = gate_by_power_hr(x, FS, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR)
    x = x(:);
    N = numel(x);

    X = fft(x);

    k = (0:N-1).';
    f = k*(FS/N);
    f_abs = min(f, FS - f);
    bpm_full = 60*f_abs;

    P2   = abs(X).^2;
    Pref = max(P2) + eps;
    Pdb  = 10*log10(P2/Pref + eps);

    if SMOOTH_BINS > 1
        Pdbs = movmean(Pdb, SMOOTH_BINS, 'omitnan');
    else
        Pdbs = Pdb;
    end

    mHR = (bpm_full >= HR_MIN_BPM) & (bpm_full <= HR_MAX_BPM);
    if ~any(mHR)
        xatt = x;
        dbg = struct('bpm',[],'Pdb_before',[],'Pdb_after',[],'Pdb_hr_peak',NaN,'Pdb_thr',NaN);
        return;
    end

    Pdb_hr_peak = max(Pdbs(mHR), [], 'omitnan');
    thr = Pdb_hr_peak + DELTA_DB;

    W = ones(N,1);
    mOUT = ~mHR;
    W(mOUT & (Pdbs > thr)) = GATE_FLOOR;

    if SMOOTH_BINS > 1
        W = movmean(W, SMOOTH_BINS, 'omitnan');
        W = min(1, max(GATE_FLOOR, W));
        W(mHR) = 1;
    end

    Xg   = X .* W;
    xatt = real(ifft(Xg, 'symmetric'));

    kk   = 0:floor(N/2);
    ff   = kk*(FS/N);
    bpmU = 60*ff;

    P_before = abs(X(1:numel(kk))).^2;
    P_after  = abs(Xg(1:numel(kk))).^2;

    if N > 1 && numel(P_before) > 2
        P_before(2:end-1) = 2*P_before(2:end-1);
        P_after(2:end-1)  = 2*P_after(2:end-1);
    end

    dbg.bpm        = bpmU(:);
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
