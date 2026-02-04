%% ========================================================================
% GRID POLAR-RADAR (3 PARES) | RAW -> (WRAP) -> SEG -> RESAMPLE -> FLOW -> CWT -> RIDGE
%
% Pares:
%   (1) POLARH2  <-> phases.csv
%   (2) POLARH3  <-> phases_raw.csv
%   (3) POLARH10 <-> phase.csv
%
% Saídas:
%   - TOP-5 por métrica (SCORE, RMSE, MAE, CORR) de cada dataset
%   - TOP-5 geral (por SCORE) de cada dataset
%   - TOP-5 agregado (média nos 3) por SCORE/RMSE/MAE/CORR
%
% Paralelo:
%   - parfor + DataQueue (status a cada PRINT_EVERY)
%
% IMPORTANTÍSSIMO (pedido):
%   - Se GATE/BP/HARM estiver OFF, as configs desse módulo NÃO variam (não entram no grid).
%     (ou seja: gate_off => só 1 "dummy config" de gate; idem bp/harm)
% ========================================================================

clc; clear; close all;

%% ===================== PARES (EDITE AQUI) =====================
DATASETS = repmat(struct( ...
    'name', '', ...
    'CSV_PATH', '', ...
    'POLAR_PATH', '', ...
    'COL_T_MS', 1, ...
    'COL_PHASE', 2), 1, 3);

DATASETS(1).name       = 'PolarH2';
DATASETS(1).CSV_PATH   = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases.csv';
DATASETS(1).POLAR_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH2.txt';

DATASETS(2).name       = 'PolarH3';
DATASETS(2).CSV_PATH   = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv';
DATASETS(2).POLAR_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH3.txt';

DATASETS(3).name       = 'PolarH10';
DATASETS(3).CSV_PATH   = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phase.csv';
DATASETS(3).POLAR_PATH = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH10.txt';

%% ===================== SCORE (AJUSTÁVEL) =====================
% Ajuste aqui se quiser priorizar outra coisa.
W_RMSE  = 1.00;
W_MAE   = 0.35;
W_1CORR = 5.00;   % penaliza (1 - corr_clamp)

%% ===================== PARALELO / STATUS =====================
USE_PARFOR  = true;
PRINT_EVERY = 200;       % status a cada N configs (no parfor, via DataQueue)
TOPK        = 5;

%% ===================== GRID (COARSE - edite) =====================
GRID.gap_thr     = [0.5];
GRID.MIN_SEG_SEC = [2.0];

GRID.FS_TARGET = [25.0];
GRID.IMETH     = {'linear'};                 % 'linear' | 'pchip'

GRID.WRAP_ON   = [0];                      % 0=unwrap-mean | 1=wrap
GRID.FLOW      = {'GATE_THEN_BP', 'BP_THEN_GATE'};

% GATE
GRID.GATE_ON     = [1];                      % pode [0 1]
GRID.HR_MIN_BPM  = [55];
GRID.HR_MAX_BPM  = [210];
GRID.DELTA_DB    = [-3 0 3];
GRID.SMOOTH_BINS = [1 3 7];
GRID.GATE_FLOOR  = [0.0 0.1];

% BP
GRID.BP_ON    = [1];                         % pode [0 1]
GRID.BP_ORDER = [4 8];
GRID.BP_WI    = [0.7 1];
GRID.BP_WF    = [2.7];

% CWT
GRID.WI       = [0.5];
GRID.WF       = [3.0];
GRID.VOICES   = [16];
GRID.CWT_WAVE = {'amor'};

% RIDGE (area + punish + event)
GRID.BAND_W_HZ      = [0.1 0.2 0.3];
GRID.DB_FLOOR_BAND  = [-25];
GRID.RIDGE_LAMBDA   = [0 0.05 0.15 0.25];
GRID.F_JUMP_HZ      = [0.10];
GRID.BURNIN_SEC     = [0];
GRID.K_EVENT        = [30];
GRID.PCHAIN_GAMMA   = [0.8];
GRID.PCHAIN_MAX     = [12];

% HARM
GRID.HARM_ON      = [0 1];                     % pode [0 1]
GRID.HARM_REL_DB  = [-6 -2];
GRID.HARM_MIN_BPM = [40];
GRID.HARM_MAX_BPM = [220];

% smooth final (dentro do segmento)
GRID.MOVMEAN_SEC = [7];

%% ===================== MATERIALIZA GRID (com lógica OFF => não varia) =====================
CFG = materialize_grid_conditional(GRID);
nCfg = numel(CFG);
nDS  = numel(DATASETS);

fprintf('Grid configs (condicional OFF): %d | Datasets: %d\n', nCfg, nDS);

%% ===================== PRÉ-CACHE DATASETS (carrega 1x) =====================
CACHE = repmat(struct( ...
    'name','', ...
    't0',[], ...
    'x_wrap',[], ...
    'x_unwrap',[], ...
    't_polar',[], ...
    'h_polar',[]), 1, nDS);

for d = 1:nDS
    ds = DATASETS(d);
    CACHE(d).name = ds.name;

    % --- radar csv ---
    [t0, ph0] = read_radar_csv(ds.CSV_PATH, ds.COL_T_MS, ds.COL_PHASE);

    t0  = (t0 - t0(1)); % segundos
    [t0, ia] = unique(t0,'stable');
    ph0 = ph0(ia);

    ph0 = force_finite_vector(double(ph0));
    phU = unwrap(ph0);

    CACHE(d).t0       = t0(:);
    CACHE(d).x_wrap   = wrap_phase(phU);
    x_unw             = phU - mean(phU,'omitnan');
    CACHE(d).x_unwrap = force_finite_vector(x_unw);

    % --- polar ---
    [tP, hP] = read_txt_polar_flex(ds.POLAR_PATH);
    CACHE(d).t_polar = tP(:);
    CACHE(d).h_polar = hP(:);
end

%% ===================== RESULT ARRAYS =====================
RMSE  = nan(nCfg, nDS);
MAE   = nan(nCfg, nDS);
CORR  = nan(nCfg, nDS);
SCORE = nan(nCfg, nDS);
NVAL  = zeros(nCfg, nDS);

%% ===================== LOOP GRID (parfor) =====================
tStart = tic;

if USE_PARFOR && license('test','Distrib_Computing_Toolbox')
    dq = parallel.pool.DataQueue;
    progress = 0;
    afterEach(dq, @(msg) local_progress_print(msg));

    parfor c = 1:nCfg
        [RMSE(c,:), MAE(c,:), CORR(c,:), SCORE(c,:), NVAL(c,:)] = ...
            eval_one_cfg(CFG(c), CACHE, W_RMSE, W_MAE, W_1CORR);

        if mod(c, PRINT_EVERY) == 0 || c == 1
            send(dq, struct('c',c,'nCfg',nCfg,'cfg_id',CFG(c).cfg_id,'elapsed',toc(tStart)));
        end
    end
else
    for c = 1:nCfg
        if mod(c, PRINT_EVERY) == 0 || c == 1
            fprintf('[%5d/%5d] cfg_id=%d | elapsed=%.1fs\n', c, nCfg, CFG(c).cfg_id, toc(tStart));
        end
        [RMSE(c,:), MAE(c,:), CORR(c,:), SCORE(c,:), NVAL(c,:)] = ...
            eval_one_cfg(CFG(c), CACHE, W_RMSE, W_MAE, W_1CORR);
    end
end

fprintf('DONE. elapsed=%.1fs\n', toc(tStart));

%% ===================== AGREGADO (média nos 3 datasets) =====================
RMSEg  = mean(RMSE,  2, 'omitnan');
MAEg   = mean(MAE,   2, 'omitnan');
CORRg  = mean(CORR,  2, 'omitnan');
SCOREg = mean(SCORE, 2, 'omitnan');
TOPK = 50;  % use 20/50/100 pra ficar estável

deltas = arrayfun(@(c) c.DELTA_DB, CFG).';  % coluna com DELTA_DB por cfg

% ----------- agregado -----------
report_dom('AGG SCOREg', deltas, SCOREg, 'ascend', TOPK);
report_dom('AGG RMSEg',  deltas, RMSEg,  'ascend', TOPK);
report_dom('AGG MAEg',   deltas, MAEg,   'ascend', TOPK);
report_dom('AGG CORRg',  deltas, CORRg,  'descend', TOPK);

% ----------- por dataset -----------
for d=1:size(SCORE,2)
    report_dom(sprintf('%s SCORE',  CACHE(d).name), deltas, SCORE(:,d), 'ascend', TOPK);
    report_dom(sprintf('%s RMSE',   CACHE(d).name), deltas, RMSE(:,d),  'ascend', TOPK);
    report_dom(sprintf('%s MAE',    CACHE(d).name), deltas, MAE(:,d),   'ascend', TOPK);
    report_dom(sprintf('%s CORR',   CACHE(d).name), deltas, CORR(:,d),  'descend', TOPK);
end

function report_dom(tag, deltas, metric, mode, TOPK)
    if strcmpi(mode,'ascend')
        [~, ix] = sort(metric, 'ascend', 'MissingPlacement','last');
    else
        [~, ix] = sort(metric, 'descend', 'MissingPlacement','last');
    end
    ix = ix(1:min(TOPK, numel(ix)));
    v  = deltas(ix);
    v  = v(isfinite(v));  % só configs com gate_on=1 (as outras viram NaN)

    if isempty(v)
        fprintf('%s: sem dados (gate_off?)\n', tag);
        return;
    end

    [u,~,g] = unique(v);
    cnt = accumarray(g,1);
    [cnt, ord] = sort(cnt,'descend');
    u = u(ord);

    fprintf('%s | TOP-%d: ', tag, TOPK);
    for i=1:numel(u)
        fprintf('DELTA=%g (%d)  ', u(i), cnt(i));
    end
    fprintf('\n');
end


%% ===================== TOP-5 AGREGADO =====================
fprintf('\n================== TOP-%d AGREGADO por SCOREg (menor melhor) ==================\n', TOPK);
ix = sort_idx(SCOREg,'ascend');
print_top_aggr(CFG, ix(1:TOPK), SCOREg, RMSEg, MAEg, CORRg);

fprintf('\n================== TOP-%d AGREGADO por RMSEg (menor melhor) ==================\n', TOPK);
ix = sort_idx(RMSEg,'ascend');
print_top_aggr(CFG, ix(1:TOPK), SCOREg, RMSEg, MAEg, CORRg);

fprintf('\n================== TOP-%d AGREGADO por MAEg (menor melhor) ==================\n', TOPK);
ix = sort_idx(MAEg,'ascend');
print_top_aggr(CFG, ix(1:TOPK), SCOREg, RMSEg, MAEg, CORRg);

fprintf('\n================== TOP-%d AGREGADO por CORRg (maior melhor) ==================\n', TOPK);
ix = sort_idx(CORRg,'descend');
print_top_aggr(CFG, ix(1:TOPK), SCOREg, RMSEg, MAEg, CORRg);

%% ===================== TOP-5 POR DATASET =====================
for d = 1:nDS
    nm = CACHE(d).name;

    fprintf('\n================== %s | TOP-%d GERAL (por SCORE) ==================\n', nm, TOPK);
    ix = sort_idx(SCORE(:,d),'ascend');
    print_top_1ds(CFG, ix(1:TOPK), SCORE(:,d), RMSE(:,d), MAE(:,d), CORR(:,d), NVAL(:,d));

    fprintf('\n================== %s | TOP-%d por RMSE ==================\n', nm, TOPK);
    ix = sort_idx(RMSE(:,d),'ascend');
    print_top_1ds(CFG, ix(1:TOPK), SCORE(:,d), RMSE(:,d), MAE(:,d), CORR(:,d), NVAL(:,d));

    fprintf('\n================== %s | TOP-%d por MAE ==================\n', nm, TOPK);
    ix = sort_idx(MAE(:,d),'ascend');
    print_top_1ds(CFG, ix(1:TOPK), SCORE(:,d), RMSE(:,d), MAE(:,d), CORR(:,d), NVAL(:,d));

    fprintf('\n================== %s | TOP-%d por CORR ==================\n', nm, TOPK);
    ix = sort_idx(CORR(:,d),'descend');
    print_top_1ds(CFG, ix(1:TOPK), SCORE(:,d), RMSE(:,d), MAE(:,d), CORR(:,d), NVAL(:,d));
end

%% ===================== (OPCIONAL) PLOT do melhor cfg agregado =====================
best = find(SCOREg == min(SCOREg,[],'omitnan'), 1, 'first');
if ~isempty(best)
    fprintf('\nPlot do melhor AGREGADO: cfg_id=%d\n', CFG(best).cfg_id);
    plot_best_cfg(CFG(best), CACHE);
end

%% ========================================================================
%                               FUNÇÕES
% ========================================================================

function local_progress_print(msg)
    fprintf('[%5d/%5d] cfg_id=%d | elapsed=%.1fs\n', msg.c, msg.nCfg, msg.cfg_id, msg.elapsed);
end

function ix = sort_idx(v, mode)
    if strcmpi(mode,'ascend')
        [~, ix] = sort(v, 'ascend', 'MissingPlacement','last');
    else
        [~, ix] = sort(v, 'descend', 'MissingPlacement','last');
    end
end

function CFG = materialize_grid_conditional(G)
    % Grid condicional:
    % - escolhe GATE_ON; se 0 => gate params viram 1 combo só (dummy)
    % - escolhe BP_ON;   se 0 => bp params viram 1 combo só (dummy)
    % - escolhe HARM_ON; se 0 => harm params viram 1 combo só (dummy)

    cfgs = [];
    id = 0;

    for gap_thr = G.gap_thr
    for minseg  = G.MIN_SEG_SEC
    for fs      = G.FS_TARGET
    for im_i    = 1:numel(G.IMETH)
        IMETH = G.IMETH{im_i};
    for wrap_on = G.WRAP_ON
    for flow_i  = 1:numel(G.FLOW)
        FLOW = G.FLOW{flow_i};

        for gate_on = G.GATE_ON
            if gate_on == 1
                HR_MIN_BPM_list  = G.HR_MIN_BPM;
                HR_MAX_BPM_list  = G.HR_MAX_BPM;
                DELTA_DB_list    = G.DELTA_DB;
                SMOOTH_BINS_list = G.SMOOTH_BINS;
                GATE_FLOOR_list  = G.GATE_FLOOR;
            else
                HR_MIN_BPM_list  = NaN;
                HR_MAX_BPM_list  = NaN;
                DELTA_DB_list    = NaN;
                SMOOTH_BINS_list = NaN;
                GATE_FLOOR_list  = NaN;
            end

            for bp_on = G.BP_ON
                if bp_on == 1
                    BP_ORDER_list = G.BP_ORDER;
                    BP_WI_list    = G.BP_WI;
                    BP_WF_list    = G.BP_WF;
                else
                    BP_ORDER_list = NaN;
                    BP_WI_list    = NaN;
                    BP_WF_list    = NaN;
                end

                for harm_on = G.HARM_ON
                    if harm_on == 1
                        HARM_REL_DB_list  = G.HARM_REL_DB;
                        HARM_MIN_BPM_list = G.HARM_MIN_BPM;
                        HARM_MAX_BPM_list = G.HARM_MAX_BPM;
                    else
                        HARM_REL_DB_list  = NaN;
                        HARM_MIN_BPM_list = NaN;
                        HARM_MAX_BPM_list = NaN;
                    end

                    for HR_MIN_BPM = HR_MIN_BPM_list
                    for HR_MAX_BPM = HR_MAX_BPM_list
                    for DELTA_DB   = DELTA_DB_list
                    for SMOOTH_BINS= SMOOTH_BINS_list
                    for GATE_FLOOR = GATE_FLOOR_list

                    for BP_ORDER = BP_ORDER_list
                    for BP_WI    = BP_WI_list
                    for BP_WF    = BP_WF_list

                    for WI = G.WI
                    for WF = G.WF
                    for VOICES = G.VOICES
                    for cw_i = 1:numel(G.CWT_WAVE)
                        CWT_WAVE = G.CWT_WAVE{cw_i};

                        for BAND_W_HZ = G.BAND_W_HZ
                        for DB_FLOOR_BAND = G.DB_FLOOR_BAND
                        for RIDGE_LAMBDA = G.RIDGE_LAMBDA
                        for F_JUMP_HZ = G.F_JUMP_HZ
                        for BURNIN_SEC = G.BURNIN_SEC
                        for K_EVENT = G.K_EVENT
                        for GAMMA = G.PCHAIN_GAMMA
                        for KMAX = G.PCHAIN_MAX
                        for HARM_REL_DB = HARM_REL_DB_list
                        for HARM_MIN_BPM = HARM_MIN_BPM_list
                        for HARM_MAX_BPM = HARM_MAX_BPM_list
                        for MOVMEAN_SEC = G.MOVMEAN_SEC

                            % sanity (evita configs inválidas quando BP_ON=1)
                            if bp_on==1
                                if ~(isfinite(BP_WI) && isfinite(BP_WF) && BP_WI>0 && BP_WF>BP_WI && BP_WF < fs/2)
                                    continue;
                                end
                            end
                            if ~(WI>0 && WF>WI && WF < fs/2)
                                continue;
                            end

                            id = id + 1;
                            c = struct();
                            c.cfg_id = id;

                            c.gap_thr     = gap_thr;
                            c.MIN_SEG_SEC = minseg;

                            c.FS_TARGET = fs;
                            c.IMETH     = IMETH;

                            c.WRAP_ON   = wrap_on;
                            c.FLOW      = FLOW;

                            c.GATE_ON     = gate_on;
                            c.HR_MIN_BPM  = HR_MIN_BPM;
                            c.HR_MAX_BPM  = HR_MAX_BPM;
                            c.DELTA_DB    = DELTA_DB;
                            c.SMOOTH_BINS = SMOOTH_BINS;
                            c.GATE_FLOOR  = GATE_FLOOR;

                            c.BP_ON    = bp_on;
                            c.BP_ORDER = BP_ORDER;
                            c.BP_WI    = BP_WI;
                            c.BP_WF    = BP_WF;

                            c.WI       = WI;
                            c.WF       = WF;
                            c.VOICES   = VOICES;
                            c.CWT_WAVE = CWT_WAVE;

                            c.BAND_W_HZ      = BAND_W_HZ;
                            c.DB_FLOOR_BAND  = DB_FLOOR_BAND;
                            c.RIDGE_LAMBDA   = RIDGE_LAMBDA;
                            c.F_JUMP_HZ      = F_JUMP_HZ;
                            c.BURNIN_SEC     = BURNIN_SEC;
                            c.K_EVENT        = K_EVENT;
                            c.PCHAIN_GAMMA   = GAMMA;
                            c.PCHAIN_MAX     = KMAX;

                            c.HARM_ON      = harm_on;
                            c.HARM_REL_DB  = HARM_REL_DB;
                            c.HARM_MIN_BPM = HARM_MIN_BPM;
                            c.HARM_MAX_BPM = HARM_MAX_BPM;

                            c.MOVMEAN_SEC = MOVMEAN_SEC;

                            cfgs = [cfgs; c]; %#ok<AGROW>

                        end; end; end; end; end; end; end; end; end; end; end; end; end; end
                    end; end; end; end
                    end; end; end; end; end
                    end; end; end
                    end; end; end; end
                end
            end
        end


    CFG = cfgs;
end

function [RMSEd, MAEd, CORRd, SCOREd, NVALd] = eval_one_cfg(cfg, CACHE, W_RMSE, W_MAE, W_1CORR)
    nDS = numel(CACHE);

    RMSEd  = nan(1,nDS);
    MAEd   = nan(1,nDS);
    CORRd  = nan(1,nDS);
    SCOREd = nan(1,nDS);
    NVALd  = zeros(1,nDS);

    for d = 1:nDS
        t0 = CACHE(d).t0;
        if isempty(t0), continue; end

        segments = segment_by_gaps(t0, cfg.gap_thr);

        if cfg.WRAP_ON == 1
            x0 = CACHE(d).x_wrap;
        else
            x0 = CACHE(d).x_unwrap;
        end

        [t_all, h_all] = radar_pipeline_segments(t0, x0, segments, cfg);

        tP = CACHE(d).t_polar;
        hP = CACHE(d).h_polar;

        if isempty(tP) || isempty(hP) || isempty(t_all)
            continue;
        end

        hC_onP = interp1(t_all, h_all, tP, 'linear', NaN);

        m = isfinite(hC_onP) & isfinite(hP);
        NVALd(d) = nnz(m);

        if nnz(m) < 5
            continue;
        end

        e = hC_onP(m) - hP(m);

        RMSEd(d) = sqrt(mean(e.^2,'omitnan'));
        MAEd(d)  = mean(abs(e),'omitnan');
        CORRd(d) = corr(hC_onP(m), hP(m), 'Rows','complete');

        corr_clamp = max(0, min(1, CORRd(d)));
        SCOREd(d) = W_RMSE*RMSEd(d) + W_MAE*MAEd(d) + W_1CORR*(1 - corr_clamp);
    end
end

function [t_all, h_all] = radar_pipeline_segments(t0, x0, segments, cfg)
    FS = cfg.FS_TARGET;
    dt = 1/FS;

    tseg_list  = {};
    hrseg_list = {};

    for s = 1:size(segments,1)
        i0 = segments(s,1); i1 = segments(s,2);
        ts = t0(i0:i1);
        xs = x0(i0:i1);

        [ts, iu] = unique(ts,'stable');
        xs = xs(iu);

        if numel(ts) < 8, continue; end

        tnew = (ts(1):dt:ts(end))';
        xnew = interp1(ts, xs, tnew, cfg.IMETH);

        ok = isfinite(tnew) & isfinite(xnew);
        tnew = tnew(ok);
        xnew = xnew(ok);

        if numel(tnew) < 16, continue; end
        if (tnew(end) - tnew(1)) < cfg.MIN_SEG_SEC, continue; end

        x_fin = force_finite_vector(xnew);

        % ================= FLOW (com GATE/BP condicional) =================
        if cfg.GATE_ON ~= 1
            % GATE OFF => só BP se ON
            if cfg.BP_ON == 1
                x_fin = apply_bp_sos(x_fin, FS, cfg.BP_ORDER, cfg.BP_WI, cfg.BP_WF);
            end
        else
            if strcmpi(cfg.FLOW,'BP_THEN_GATE')
                if cfg.BP_ON == 1
                    x_fin = apply_bp_sos(x_fin, FS, cfg.BP_ORDER, cfg.BP_WI, cfg.BP_WF);
                end
                x_fin = gate_by_power_hr_full(x_fin, FS, cfg.HR_MIN_BPM, cfg.HR_MAX_BPM, cfg.DELTA_DB, cfg.SMOOTH_BINS, cfg.GATE_FLOOR);
            else
                x_fin = gate_by_power_hr_full(x_fin, FS, cfg.HR_MIN_BPM, cfg.HR_MAX_BPM, cfg.DELTA_DB, cfg.SMOOTH_BINS, cfg.GATE_FLOOR);
                if cfg.BP_ON == 1
                    x_fin = apply_bp_sos(x_fin, FS, cfg.BP_ORDER, cfg.BP_WI, cfg.BP_WF);
                end
            end
        end

        x_fin = force_finite_vector(x_fin);

        % ================= CWT =================
        [wt, fbin] = cwt(x_fin, cfg.CWT_WAVE, FS, ...
            'FrequencyLimits',[cfg.WI cfg.WF], ...
            'VoicesPerOctave', cfg.VOICES);

        fbin = fbin(:);
        if numel(fbin)>1 && fbin(2) < fbin(1)
            fbin = flipud(fbin);
            wt   = flipud(wt);
        end

        P = single(abs(wt).^2);
        clear wt;
        P(~isfinite(P)) = 0;

        % ================= RIDGE (com HARM condicional) =================
        freq_hz = ridge_area_punish_event_harm(P, fbin, FS, ...
            cfg.BAND_W_HZ, cfg.DB_FLOOR_BAND, cfg.RIDGE_LAMBDA, cfg.F_JUMP_HZ, cfg.BURNIN_SEC, cfg.K_EVENT, cfg.PCHAIN_GAMMA, cfg.PCHAIN_MAX, ...
            cfg.HARM_ON, cfg.HARM_REL_DB, cfg.HARM_MIN_BPM, cfg.HARM_MAX_BPM);

        hrseg = 60*freq_hz(:);

        if cfg.MOVMEAN_SEC > 0
            Wseg = max(1, round(cfg.MOVMEAN_SEC * FS));
            hrseg = movmean(hrseg, Wseg, 'Endpoints','shrink');
        end

        tseg_list{end+1}  = tnew(:); %#ok<AGROW>
        hrseg_list{end+1} = hrseg(:); %#ok<AGROW>
    end

    [t_all, h_all] = concat_segments(tseg_list, hrseg_list);
end

function print_top_aggr(CFG, idx, SCOREg, RMSEg, MAEg, CORRg)
    for k = 1:numel(idx)
        i = idx(k);
        c = CFG(i);
        fprintf('#%02d cfg_id=%d | SCOREg=%.4f RMSEg=%.4f MAEg=%.4f CORRg=%.4f | WRAP=%d FLOW=%s G=%d BP=%d H=%d | BP[%.2f %.2f] ord=%g | lam=%.3f jump=%.3f gamma=%.2f kmax=%g bw=%.2f\n', ...
            k, c.cfg_id, SCOREg(i), RMSEg(i), MAEg(i), CORRg(i), ...
            c.WRAP_ON, string(c.FLOW), c.GATE_ON, c.BP_ON, c.HARM_ON, ...
            c.BP_WI, c.BP_WF, c.BP_ORDER, ...
            c.RIDGE_LAMBDA, c.F_JUMP_HZ, c.PCHAIN_GAMMA, c.PCHAIN_MAX, c.BAND_W_HZ);
    end
end

function print_top_1ds(CFG, idx, SCORE, RMSE, MAE, CORR, NVAL)
    for k = 1:numel(idx)
        i = idx(k);
        c = CFG(i);
        fprintf('#%02d cfg_id=%d | SCORE=%.4f RMSE=%.4f MAE=%.4f CORR=%.4f N=%d | WRAP=%d FLOW=%s G=%d BP=%d H=%d | BP[%.2f %.2f] ord=%g | lam=%.3f jump=%.3f gamma=%.2f kmax=%g bw=%.2f\n', ...
            k, c.cfg_id, SCORE(i), RMSE(i), MAE(i), CORR(i), NVAL(i), ...
            c.WRAP_ON, string(c.FLOW), c.GATE_ON, c.BP_ON, c.HARM_ON, ...
            c.BP_WI, c.BP_WF, c.BP_ORDER, ...
            c.RIDGE_LAMBDA, c.F_JUMP_HZ, c.PCHAIN_GAMMA, c.PCHAIN_MAX, c.BAND_W_HZ);
    end
end

function plot_best_cfg(cfg, CACHE)
    figure('Color','w'); tiledlayout(numel(CACHE),1,'Padding','compact','TileSpacing','compact');
    for d=1:numel(CACHE)
        nexttile; hold on; grid on;
        t0 = CACHE(d).t0;

        segments = segment_by_gaps(t0, cfg.gap_thr);
        if cfg.WRAP_ON==1, x0=CACHE(d).x_wrap; else, x0=CACHE(d).x_unwrap; end

        [t_all, h_all] = radar_pipeline_segments(t0, x0, segments, cfg);

        tP = CACHE(d).t_polar;
        hP = CACHE(d).h_polar;

        if isempty(tP) || isempty(hP) || isempty(t_all)
            title(sprintf('%s | sem dados', CACHE(d).name)); continue;
        end

        hC_onP = interp1(t_all, h_all, tP, 'linear', NaN);

        plot(tP, hP, 'c--', 'LineWidth', 1.2);
        plot(tP, hC_onP, 'k-', 'LineWidth', 1.2);
        xlabel('t (s)'); ylabel('HR (bpm)');
        title(sprintf('%s | cfg_id=%d', CACHE(d).name, cfg.cfg_id));
        legend('POLAR','RADAR@POLAR','Location','best');
    end
end

%% ===================== LEITURA RADAR CSV =====================
function [t_sec, phase] = read_radar_csv(path, col_t_ms, col_phase)
    A = readmatrix(path);
    if isempty(A) || size(A,2) < max(col_t_ms, col_phase)
        t_sec = []; phase = []; return;
    end
    tpuro_ms = double(A(:, col_t_ms));
    phase    = double(A(:, col_phase));
    tpuro_ms = tpuro_ms(:);
    phase    = phase(:);
    ok = isfinite(tpuro_ms) & isfinite(phase);
    tpuro_ms = tpuro_ms(ok);
    phase    = phase(ok);
    t_sec = (tpuro_ms - tpuro_ms(1))/1000;
end

%% ===================== SEGMENTAÇÃO / CONCAT =====================
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
    if isempty(t_all)
        return;
    end
    [t_all, iu] = unique(t_all, 'stable');
    h_all = h_all(iu);
end

%% ===================== WRAP / SANITIZE =====================
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
        x(:) = 0;
        return;
    end

    x(~isfinite(x)) = 0;
    if ~any(isfinite(x)), x(:) = 0; end
end

%% ===================== BANDPASS (SOS) =====================
function x = apply_bp_sos(x, Fs, ord, wi, wf)
    x = force_finite_vector(x);
    if ~isfinite(Fs) || Fs<=0, return; end
    if ~isfinite(ord) || ord<1, return; end
    if ~isfinite(wi) || ~isfinite(wf), return; end
    if wf >= Fs/2, return; end
    if wi <= 0 || wi >= wf, return; end

    Wn = [wi wf]/(Fs/2);

    use_sos = false;
    try
        [sosBP, gBP] = butter(ord, Wn, 'bandpass', 'sos');
        use_sos = true;
    catch
        use_sos = false;
    end

    if ~use_sos
        [z,p,k] = butter(ord, Wn, 'bandpass');
        sosBP = zp2sos(z,p,k);
        gBP = 1;
    end

    if exist('sosfiltfilt','file') == 2
        x = sosfiltfilt(sosBP, gBP, x);
    else
        [b,a] = sos2tf(sosBP, gBP);
        x = filtfilt(b,a,x);
    end

    x = force_finite_vector(x);
end

%% ===================== GATE (por potência fora da banda HR) =====================
function xatt = gate_by_power_hr_full(x, FS, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR)
    x = force_finite_vector(x);
    N = numel(x);
    if N < 8 || ~isfinite(FS) || FS<=0
        xatt = x; return;
    end
    if ~isfinite(HR_MIN_BPM) || ~isfinite(HR_MAX_BPM) || HR_MIN_BPM>=HR_MAX_BPM
        xatt = x; return;
    end
    if ~isfinite(DELTA_DB), DELTA_DB = -3; end
    if ~isfinite(SMOOTH_BINS), SMOOTH_BINS = 7; end
    if ~isfinite(GATE_FLOOR), GATE_FLOOR = 0.0; end

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

%% ===================== RIDGE (area + punish + event + harm) =====================
function freq_hz = ridge_area_punish_event_harm(P, f_bin, srate, BAND_W_HZ, DB_FLOOR_BAND, RIDGE_LAMBDA, F_JUMP_HZ, BURNIN_SEC, K_EVENT, GAMMA, KMAX, ...
                                                HARM_ON, HARM_REL_DB, HARM_MIN_BPM, HARM_MAX_BPM)
    f_bin = double(f_bin(:));
    nb = numel(f_bin);
    nt = size(P,2);

    if nb < 2 || nt < 2
        freq_hz = nan(nt,1); return;
    end

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

%% ===================== POLAR READER =====================
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
        try
            t_dt = datetime(ts_str,'InputFormat',"yyyy-MM-dd HH:mm:ss.SSS");
        catch
            t_dt = datetime(ts_str);
        end
    end

    t_sec = seconds(t_dt - t_dt(1));
    t_sec = t_sec(:); HR = HR(:);

    m = isfinite(t_sec) & isfinite(HR);
    t_sec = t_sec(m); HR = HR(m);

    [t_sec, iu] = unique(t_sec,'stable');
    HR = HR(iu);
end
