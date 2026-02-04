clc; clear; close all;

%% ===================== PATHS (AJUSTE) =====================
BASE = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\';

DSETS = struct([]);

DSETS(1).name       = 'POLARH2 <-> phases.csv';
DSETS(1).csv_path   = fullfile(BASE,'phases.csv');
DSETS(1).polar_path = fullfile(BASE,'POLARH2.txt');
DSETS(1).col_t_ms   = 1;  DSETS(1).col_phase = 2;

DSETS(2).name       = 'POLARH3 <-> phases_raw.csv';
DSETS(2).csv_path   = fullfile(BASE,'phases_raw.csv');
DSETS(2).polar_path = fullfile(BASE,'POLARH3.txt');
DSETS(2).col_t_ms   = 1;  DSETS(2).col_phase = 2;

DSETS(3).name       = 'POLARH10 <-> phase.csv';
DSETS(3).csv_path   = fullfile(BASE,'phase.csv');
DSETS(3).polar_path = fullfile(BASE,'POLARH10.txt');
DSETS(3).col_t_ms   = 1;  DSETS(3).col_phase = 2;

%% ===================== GRID (COARSE, REAL) =====================
% ---- segmentação / resample ----
GRID.gap_thr     = [0.5];      % s
GRID.MIN_SEG_SEC = [2.0];      % s
GRID.FS_TARGET   = [25];            % Hz
GRID.IMETH       = {'linear'}; % interp
GRID.WRAP_ON     = [0 1];              % wrap+mean vs mean

% ---- ordem do fluxo ----
GRID.FLOW = {'BP_THEN_GATE','GATE_THEN_BP'};

% ---- features ON/OFF ----
GRID.BP_ON   = [1];
GRID.GATE_ON = [1];
GRID.HARM_ON = [0 1];

% ---- BP params (só usados se BP_ON=1) ----
GRID.BP_ORDER = [2 6];
GRID.BP_WI    = [0.1 0.4 0.7];  % Hz (sem clamp)
GRID.BP_WF    = [2.5];     % Hz

% ---- GATE params (só usados se GATE_ON=1) ----
GRID.HR_MIN_BPM  = [65];
GRID.HR_MAX_BPM  = [210];
GRID.DELTA_DB    = [-3 -5];
GRID.SMOOTH_BINS = [0 3 7];          % 0 = off
GRID.GATE_FLOOR  = [0 0.1];

% ---- CWT (mantém quase fixo pra não explodir) ----
GRID.WI       = [0.6];
GRID.WF       = [3.0];
GRID.VOICES   = [16];
GRID.CWT_WAVE = {'amor'};

% ---- RIDGE (varre o que costuma mandar muito) ----
GRID.BAND_W_HZ     = [0.15 0.30];
GRID.DB_FLOOR_BAND = [-25];
GRID.RIDGE_LAMBDA  = [0.05 0.15 0.3];
GRID.F_JUMP_HZ     = [0.15 0.4];
GRID.BURNIN_SEC    = [5];
GRID.K_EVENT       = [30];
GRID.PCHAIN_GAMMA  = [0.4 0.8];
GRID.PCHAIN_MAX    = [6];

% ---- HARM params (só usados se HARM_ON=1) ----
GRID.HARM_REL_DB  = [-2 -6];
GRID.HARM_MIN_BPM = [65];
GRID.HARM_MAX_BPM = [220];

% ---- suavização final ----
GRID.MOVMEAN_SEC = [7];

%% ===================== BUILD CONFIGS (SEM CONFIGS OCIOSAS) =====================
CFG = build_cfg_list(GRID);
NCFG = numel(CFG);

fprintf('Datasets=%d | Configs=%d | Total jobs=%d\n', numel(DSETS), NCFG, numel(DSETS)*NCFG);

%% ===================== PARALLEL =====================
p = gcp('nocreate');
if isempty(p)
    parpool('Processes',6); % ou 'Threads'
end

dq = parallel.pool.DataQueue;
afterEach(dq, @progress_cb);
progress_cb('init', numel(DSETS)*NCFG);

%% ===================== LOOP DATASETS =====================
ALL = cell(numel(DSETS),1);

for d = 1:numel(DSETS)
    ds = DSETS(d);
    fprintf('\n==================== %s ====================\n', ds.name);

    % --- lê radar 1x por dataset ---
    A = readmatrix(ds.csv_path);
    tpuro_ms = A(:, ds.col_t_ms);
    phase0   = A(:, ds.col_phase);

    tpuro_ms = tpuro_ms(:);
    phase0   = phase0(:);

    t0 = (tpuro_ms - tpuro_ms(1))/1000;
    [t0, ia] = unique(t0, 'stable');
    phase0   = phase0(ia);

    ph = force_finite_vector(double(phase0));
    ph = unwrap(ph);

    % --- lê polar 1x por dataset ---
    [tG, hP] = read_txt_polar_flex(ds.polar_path);
    tG = tG(:); hP = hP(:);

    R = repmat(struct( ...
        'ds', ds.name, 'cfg_id',0, ...
        'RMSE',NaN,'MAE',NaN,'CORR',NaN,'Nval',0, ...
        'bp',0,'gate',0,'harm',0, ...
        'gap_thr',NaN,'minseg',NaN,'fs',NaN,'wrap',NaN,'imeth',"", 'flow',"", ...
        'bp_ord',NaN,'bp_wi',NaN,'bp_wf',NaN, ...
        'hrmin',NaN,'hrmax',NaN,'delta_db',NaN,'smooth_bins',NaN,'gate_floor',NaN, ...
        'voices',NaN,'bandw',NaN,'dbfloor',NaN,'lam',NaN,'fjump',NaN,'burn',NaN,'kevent',NaN,'pgamma',NaN,'pmax',NaN, ...
        'hrel',NaN,'hmin',NaN,'hmax',NaN, ...
        'mov',NaN), NCFG, 1);

    parfor c = 1:NCFG
        C = CFG(c);

        % prepara x0 conforme WRAP_ON
        if C.WRAP_ON == 1
            x0 = wrap_phase(ph);
        else
            x0 = ph - mean(ph,'omitnan');
            x0 = force_finite_vector(x0);
        end

        segments = segment_by_gaps(t0, C.gap_thr);

        [t_all, h_all] = run_one_cfg( ...
            t0, x0, segments, ...
            C.FS_TARGET, char(C.IMETH), C.MIN_SEG_SEC, ...
            char(C.FLOW), ...
            C.BP_ON, C.BP_ORDER, C.BP_WI, C.BP_WF, ...
            C.GATE_ON, C.HR_MIN_BPM, C.HR_MAX_BPM, C.DELTA_DB, C.SMOOTH_BINS, C.GATE_FLOOR, ...
            C.WI, C.WF, C.VOICES, char(C.CWT_WAVE), ...
            C.BAND_W_HZ, C.DB_FLOOR_BAND, C.RIDGE_LAMBDA, C.F_JUMP_HZ, C.BURNIN_SEC, ...
            C.K_EVENT, C.PCHAIN_GAMMA, C.PCHAIN_MAX, ...
            C.HARM_ON, C.HARM_REL_DB, C.HARM_MIN_BPM, C.HARM_MAX_BPM, ...
            C.MOVMEAN_SEC);

        hC_onP = interp1(t_all, h_all, tG, 'linear', NaN);
        mE = isfinite(hC_onP) & isfinite(hP);
        Nval = nnz(mE);

        if Nval >= 5
            e = hC_onP(mE) - hP(mE);
            RMSE = sqrt(mean(e.^2,'omitnan'));
            MAE  = mean(abs(e),'omitnan');
            CORR = corr(hC_onP(mE), hP(mE), 'Rows','complete');
        else
            RMSE = NaN; MAE = NaN; CORR = NaN;
        end

tmp = R(c);  % <<< usa template pré-alocado (mesmos campos sempre)

tmp.ds     = ds.name;
tmp.cfg_id = c;

tmp.RMSE = RMSE;
tmp.MAE  = MAE;
tmp.CORR = CORR;
tmp.Nval = Nval;

% flags
tmp.bp   = C.BP_ON;
tmp.gate = C.GATE_ON;
tmp.harm = C.HARM_ON;

% cfg dump (tudo numérico/char; evita "string" no parfor)
tmp.gap_thr = C.gap_thr;
tmp.minseg  = C.MIN_SEG_SEC;
tmp.fs      = C.FS_TARGET;
tmp.wrap    = C.WRAP_ON;
tmp.imeth   = char(C.IMETH);   % <<< char
tmp.flow    = char(C.FLOW);    % <<< char

tmp.bp_ord = C.BP_ORDER;
tmp.bp_wi  = C.BP_WI;
tmp.bp_wf  = C.BP_WF;

tmp.hrmin      = C.HR_MIN_BPM;
tmp.hrmax      = C.HR_MAX_BPM;
tmp.delta_db   = C.DELTA_DB;
tmp.smooth_bins= C.SMOOTH_BINS;
tmp.gate_floor = C.GATE_FLOOR;

tmp.voices = C.VOICES;

tmp.bandw  = C.BAND_W_HZ;
tmp.dbfloor= C.DB_FLOOR_BAND;
tmp.lam    = C.RIDGE_LAMBDA;
tmp.fjump  = C.F_JUMP_HZ;
tmp.burn   = C.BURNIN_SEC;
tmp.kevent = C.K_EVENT;
tmp.pgamma = C.PCHAIN_GAMMA;
tmp.pmax   = C.PCHAIN_MAX;

tmp.hrel = C.HARM_REL_DB;
tmp.hmin = C.HARM_MIN_BPM;
tmp.hmax = C.HARM_MAX_BPM;

tmp.mov = C.MOVMEAN_SEC;

R(c) = tmp;  % <<< agora sim: struct "fatiável" pro parfor

        send(dq, 1);
    end

    ALL{d} = R;

    T = struct2table(R);
    T = T(isfinite(T.MAE),:);
    T_mae  = sortrows(T, {'MAE','RMSE'}, {'ascend','ascend'});
    T_corr = sortrows(T, {'CORR','MAE'}, {'descend','ascend'});

    fprintf('TOP-5 por MAE:\n');
    disp(T_mae(1:min(5,height(T_mae)), {'cfg_id','bp','gate','harm','MAE','RMSE','CORR','Nval','gap_thr','minseg','fs','wrap','imeth','flow','bp_ord','bp_wi','bp_wf','delta_db','gate_floor','lam','fjump','bandw','kevent','mov'}));

    fprintf('TOP-5 por CORR:\n');
    disp(T_corr(1:min(5,height(T_corr)), {'cfg_id','bp','gate','harm','MAE','RMSE','CORR','Nval','gap_thr','minseg','fs','wrap','imeth','flow','bp_ord','bp_wi','bp_wf','delta_db','gate_floor','lam','fjump','bandw','kevent','mov'}));
end

%% ===================== AGREGADO ENTRE DATASETS =====================
AG = [];
for d = 1:numel(ALL)
    AG = [AG; struct2table(ALL{d})]; %#ok<AGROW>
end

% agrupa por chaves principais do cfg
keys = {'bp','gate','harm','gap_thr','minseg','fs','wrap','imeth','flow', ...
        'bp_ord','bp_wi','bp_wf','hrmin','hrmax','delta_db','smooth_bins','gate_floor', ...
        'voices','bandw','dbfloor','lam','fjump','burn','kevent','pgamma','pmax', ...
        'hrel','hmin','hmax','mov'};

G = groupsummary(AG, keys, 'mean', {'MAE','RMSE','CORR','Nval'});
G = sortrows(G, {'mean_MAE','mean_RMSE'}, {'ascend','ascend'});

fprintf('\n==================== AGREGADO (média) | TOP-15 por MAE ====================\n');
disp(G(1:min(15,height(G)), [keys, {'mean_MAE','mean_RMSE','mean_CORR','mean_Nval'}]));

%% ===================== CALLBACK PROGRESS =====================
function progress_cb(x, total)
    persistent cnt T0 TOT
    if ischar(x) && strcmpi(x,'init')
        cnt = 0; T0 = tic; TOT = total; return;
    end
    cnt = cnt + 1;
    if mod(cnt, 200) == 0 || cnt == TOT
        fprintf('[%d/%d] elapsed=%.1fs\n', cnt, TOT, toc(T0));
    end
end

%% ===================== BUILD CFG LIST (SEM OCIOSAS) =====================
function CFG = build_cfg_list(GRID)
    CFG = struct([]);
    id = 0;

    for gap_thr = GRID.gap_thr
    for MIN_SEG_SEC = GRID.MIN_SEG_SEC
    for FS_TARGET = GRID.FS_TARGET
    for im = 1:numel(GRID.IMETH); IMETH = GRID.IMETH{im};
    for WRAP_ON = GRID.WRAP_ON
    for fl = 1:numel(GRID.FLOW); FLOW = GRID.FLOW{fl};

        for BP_ON = GRID.BP_ON
        for GATE_ON = GRID.GATE_ON
        for HARM_ON = GRID.HARM_ON

            % --- parâmetros que dependem de feature ON ---
            if BP_ON == 1
                BP_ORDER_list = GRID.BP_ORDER;
                BP_WI_list    = GRID.BP_WI;
                BP_WF_list    = GRID.BP_WF;
            else
                BP_ORDER_list = 2;  % dummy fixo
                BP_WI_list    = 0.1;
                BP_WF_list    = 3.0;
            end

            if GATE_ON == 1
                HR_MIN_list   = GRID.HR_MIN_BPM;
                HR_MAX_list   = GRID.HR_MAX_BPM;
                DELTA_list    = GRID.DELTA_DB;
                SMOOTH_list   = GRID.SMOOTH_BINS;
                FLOOR_list    = GRID.GATE_FLOOR;
            else
                HR_MIN_list = 65; HR_MAX_list = 210; DELTA_list = -5; SMOOTH_list = 0; FLOOR_list = 0.05;
            end

            if HARM_ON == 1
                HREL_list  = GRID.HARM_REL_DB;
                HMIN_list  = GRID.HARM_MIN_BPM;
                HMAX_list  = GRID.HARM_MAX_BPM;
            else
                HREL_list = -4; HMIN_list = 65; HMAX_list = 220;
            end

            for BP_ORDER = BP_ORDER_list
            for BP_WI = BP_WI_list
            for BP_WF = BP_WF_list
                if BP_ON==1 && ~(BP_WI > 0 && BP_WI < BP_WF && BP_WF < FS_TARGET/2)
                    continue;
                end

                for HR_MIN_BPM = HR_MIN_list
                for HR_MAX_BPM = HR_MAX_list
                    if GATE_ON==1 && ~(HR_MIN_BPM < HR_MAX_BPM)
                        continue;
                    end

                    for DELTA_DB = DELTA_list
                    for SMOOTH_BINS = SMOOTH_list
                    for GATE_FLOOR  = FLOOR_list

                        for WI = GRID.WI
                        for WF = GRID.WF
                            if ~(WI > 0 && WI < WF && WF < FS_TARGET/2)
                                continue;
                            end

                            for VOICES = GRID.VOICES
                            for cw = 1:numel(GRID.CWT_WAVE); CWT_WAVE = GRID.CWT_WAVE{cw};

                                for BAND_W_HZ = GRID.BAND_W_HZ
                                for DB_FLOOR_BAND = GRID.DB_FLOOR_BAND
                                for RIDGE_LAMBDA = GRID.RIDGE_LAMBDA
                                for F_JUMP_HZ = GRID.F_JUMP_HZ
                                for BURNIN_SEC = GRID.BURNIN_SEC
                                for K_EVENT = GRID.K_EVENT
                                for PCHAIN_GAMMA = GRID.PCHAIN_GAMMA
                                for PCHAIN_MAX = GRID.PCHAIN_MAX
                                for HARM_REL_DB = HREL_list
                                for HARM_MIN_BPM = HMIN_list
                                for HARM_MAX_BPM = HMAX_list
                                for MOVMEAN_SEC = GRID.MOVMEAN_SEC

                                    id = id + 1;
                                    C = struct();
                                    C.cfg_id = id;

                                    C.gap_thr = gap_thr;
                                    C.MIN_SEG_SEC = MIN_SEG_SEC;
                                    C.FS_TARGET = FS_TARGET;
                                    C.IMETH = IMETH;
                                    C.WRAP_ON = WRAP_ON;
                                    C.FLOW = FLOW;

                                    C.BP_ON = BP_ON;
                                    C.GATE_ON = GATE_ON;
                                    C.HARM_ON = HARM_ON;

                                    C.BP_ORDER = BP_ORDER;
                                    C.BP_WI = BP_WI;
                                    C.BP_WF = BP_WF;

                                    C.HR_MIN_BPM  = HR_MIN_BPM;
                                    C.HR_MAX_BPM  = HR_MAX_BPM;
                                    C.DELTA_DB    = DELTA_DB;
                                    C.SMOOTH_BINS = SMOOTH_BINS;
                                    C.GATE_FLOOR  = GATE_FLOOR;

                                    C.WI = WI; C.WF = WF;
                                    C.VOICES = VOICES;
                                    C.CWT_WAVE = CWT_WAVE;

                                    C.BAND_W_HZ = BAND_W_HZ;
                                    C.DB_FLOOR_BAND = DB_FLOOR_BAND;
                                    C.RIDGE_LAMBDA = RIDGE_LAMBDA;
                                    C.F_JUMP_HZ = F_JUMP_HZ;
                                    C.BURNIN_SEC = BURNIN_SEC;
                                    C.K_EVENT = K_EVENT;
                                    C.PCHAIN_GAMMA = PCHAIN_GAMMA;
                                    C.PCHAIN_MAX = PCHAIN_MAX;

                                    C.HARM_REL_DB = HARM_REL_DB;
                                    C.HARM_MIN_BPM = HARM_MIN_BPM;
                                    C.HARM_MAX_BPM = HARM_MAX_BPM;

                                    C.MOVMEAN_SEC = MOVMEAN_SEC;

                                    CFG = [CFG; C]; %#ok<AGROW>
                                end
                                end
                                end
                                end
                                end
                                end
                                end
                                end

                            end
                            end
                        end
                        end

                    end
                    end
                    end
                    end
                end
                end

            end
            end
            end

        end
        end
        end

    end
    end
    end
    end
    end
    end

    fprintf('CFG built: %d configs\n', numel(CFG));
end
    end
    end
end
%% ===================== RUN ONE CONFIG (PIPELINE) =====================
function [t_all, h_all] = run_one_cfg( ...
    t0, x0, segments, ...
    FS_TARGET, IMETH, MIN_SEG_SEC, ...
    FLOW, ...
    BP_ON, BP_ORDER, BP_WI, BP_WF, ...
    GATE_ON, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR, ...
    WI, WF, VOICES, CWT_WAVE, ...
    BAND_W_HZ, DB_FLOOR_BAND, RIDGE_LAMBDA, F_JUMP_HZ, BURNIN_SEC, ...
    K_EVENT, PCHAIN_GAMMA, PCHAIN_MAX, ...
    HARM_ON, HARM_REL_DB, HARM_MIN_BPM, HARM_MAX_BPM, ...
    MOVMEAN_SEC)

    dt = 1/FS_TARGET;

    tseg_list  = {};
    hrseg_list = {};

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

        % ---------- FLOW (BP/GATE) com dummies ----------
        x_fin = xnew;

        if (GATE_ON ~= 1) && (BP_ON ~= 1)
            % dummy total
        elseif GATE_ON ~= 1
            x_fin = apply_bp_sos_unitgain(x_fin, FS_TARGET, BP_ORDER, BP_WI, BP_WF, BP_ON);
        elseif BP_ON ~= 1
            x_fin = gate_by_power_hr_full(x_fin, FS_TARGET, ...
                HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR, GATE_ON);
        else
            if strcmpi(FLOW,'BP_THEN_GATE')
                x_fin = apply_bp_sos_unitgain(x_fin, FS_TARGET, BP_ORDER, BP_WI, BP_WF, BP_ON);
                x_fin = gate_by_power_hr_full(x_fin, FS_TARGET, ...
                    HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR, GATE_ON);
            else
                x_fin = gate_by_power_hr_full(x_fin, FS_TARGET, ...
                    HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR, GATE_ON);
                x_fin = apply_bp_sos_unitgain(x_fin, FS_TARGET, BP_ORDER, BP_WI, BP_WF, BP_ON);
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

        % ---------- RIDGE (+HARM dummy) ----------
        freq_hz = ridge_area_punish_event_harm(P, fbin, FS_TARGET, ...
            BAND_W_HZ, DB_FLOOR_BAND, RIDGE_LAMBDA, F_JUMP_HZ, BURNIN_SEC, ...
            K_EVENT, PCHAIN_GAMMA, PCHAIN_MAX, ...
            HARM_ON, HARM_REL_DB, HARM_MIN_BPM, HARM_MAX_BPM);

        hrseg = 60*freq_hz(:);

        if MOVMEAN_SEC > 0
            Wseg = max(1, round(MOVMEAN_SEC * FS_TARGET));
            hrseg = movmean(hrseg, Wseg, 'Endpoints','shrink');
        end

        tseg_list{end+1}  = tnew(:); %#ok<AGROW>
        hrseg_list{end+1} = hrseg(:); %#ok<AGROW>
    end

    [t_all, h_all] = concat_segments(tseg_list, hrseg_list);
end

%% ===================== FUNÇÕES BASE (iguais/compatíveis com seu code) =====================

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

% ======= BP com ganho unitário (DUMMY quando BP_ON=0) =======
function x = apply_bp_sos_unitgain(x, Fs, ord, wi, wf, BP_ON)
    x = force_finite_vector(x);

    if BP_ON ~= 1
        return; % dummy
    end

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

% ======= GATE (DUMMY quando GATE_ON=0) =======
function xatt = gate_by_power_hr_full(x, FS, HR_MIN_BPM, HR_MAX_BPM, DELTA_DB, SMOOTH_BINS, GATE_FLOOR, GATE_ON)
    x = force_finite_vector(x);

    if GATE_ON ~= 1
        xatt = x; return; % dummy
    end

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

% ===== ridge (o seu) =====
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
