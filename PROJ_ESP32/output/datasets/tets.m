clc; clear; close all;

%% ===================== BASE DIR (AJUSTE AQUI) ===================== %%
BASE_DIR = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\datasets\Figshare\dta';

% Pastas alvo (pula se não existir)
FOLDER_LIST = arrayfun(@(k) sprintf('GDN%04d',k), [1 2 3], 'UniformOutput', false);

% Padrão de arquivos por pasta
FILE_GLOB_FMT = '%s_*.mat';

%% ===================== STATUS / CONTROLE ===================== %%
VERBOSE_STATUS  = 1;
PAR_STATUS      = 1;
PAR_PRINT_EVERY = 20;

%% ===================== PARALLEL (PROCESSES) ===================== %%
MAX_WORKERS_SAFE = 12;   % ajuste conforme RAM

%% ===================== GRID (FOCO: GATE / BP / HARM + ORDEM) ===================== %%
GRID_INPUT_MODE   = {'radar'};
GRID_SRATE_TARGET = [25];

% CWT / banda de busca (fixo por enquanto)
GRID_WI       = [0.6];
GRID_WF       = [3];
GRID_VOICES   = [16];
GRID_CWT_WAVE = {'amor'};

% --------- GATE ----------
GRID_GATE_ON     = [1 0];     % ativa/desativa
GRID_HR_MIN_BPM  = [55 65];
GRID_HR_MAX_BPM  = [210];
GRID_DELTA_DB    = [-6 -3 0]; % se GATE_ON=0 será ignorado (poupado)
GRID_SMOOTH_BINS = [1 3 7];   % idem
GRID_GATE_FLOOR  = [0.0 0.1];% idem

% --------- BP ----------
GRID_BP_ON    = [1];          % se quiser testar sem BP, adicione 0 aqui
GRID_BP_ORDER = [4 8];
GRID_BP_WI    = [0.4 0.7 1.0];
GRID_BP_WF    = [2.5];

% --------- ORDEM (comparação BP antes/depois do gate) ----------
% Só faz sentido quando GATE_ON=1 e BP_ON=1.
GRID_FLOW = {'GATE_THEN_BP','BP_THEN_GATE'};

% --------- RIDGE / EVENT / PUNISH (fixo; ajuste se quiser) ----------
GRID_BAND_W_HZ       = [0.15];
GRID_DB_FLOOR_BAND   = [-25];
GRID_RIDGE_LAMBDA    = [0.25];
GRID_F_JUMP_HZ       = [0.15];
GRID_BURNIN_SEC      = [0];
GRID_K_EVENT         = [30];
GRID_PCHAIN_GAMMA    = [0.8];
GRID_PCHAIN_MAX      = [12];

% --------- HARM (f/2) ----------
GRID_HARM_ON      = [1 0];       % ativa/desativa
GRID_HARM_REL_DB  = [-8 -4]; % se HARM_ON=0 será ignorado (poupado)
GRID_HARM_MIN_BPM = [40];        % idem
GRID_HARM_MAX_BPM = [220];       % idem

% --------- ECG / avaliação (fixo) ----------
GRID_ECG_MIN_DIST_SEC = [0.35];
GRID_ECG_MIN_PROM     = [2.0];
GRID_HR_MIN_BPM_ECG   = [30];
GRID_HR_MAX_BPM_ECG   = [220];
GRID_WIN_MEAN_SEC     = [7];

%% ===================== ENUMERA ARQUIVOS ===================== %%
files = {};
ds_id  = {};
for i = 1:numel(FOLDER_LIST)
    ds = FOLDER_LIST{i};
    pds = fullfile(BASE_DIR, ds);
    if ~exist(pds,'dir'), continue; end
    D = dir(fullfile(pds, sprintf(FILE_GLOB_FMT, ds)));
    for k = 1:numel(D)
        files{end+1,1} = fullfile(D(k).folder, D(k).name); %#ok<AGROW>
        ds_id{end+1,1} = ds; %#ok<AGROW>
    end
end

fprintf('Arquivos encontrados: %d\n', numel(files));
if isempty(files), error('Nenhum arquivo encontrado. Verifique BASE_DIR / padrão.'); end

% -------- identifica tipo pelo nome: _1_, _2_, _3_ etc ----------
file_type_num = nan(numel(files),1);
for i = 1:numel(files)
    [~,nm,~] = fileparts(files{i});
    file_type_num(i) = parse_type_num(nm);
end
file_type_num(~isfinite(file_type_num)) = 0;

uniqTypes = unique(file_type_num(:).', 'stable');
nTypes = numel(uniqTypes);

typeTag = strings(1,nTypes);
for t = 1:nTypes
    if uniqTypes(t)==0, typeTag(t) = "tUNK";
    else,              typeTag(t) = "t" + string(uniqTypes(t));
    end
end

type_idx = zeros(numel(files),1);
for i = 1:numel(files)
    type_idx(i) = find(uniqTypes == file_type_num(i), 1, 'first');
end

%% ===================== PRÉ-LOAD ===================== %%
raw = cell(numel(files),1);
for i = 1:numel(files)
    raw{i} = load_one_mat(files{i});
end

%% ===================== PRÉ-CACHE POR SRATE_TARGET ===================== %%
uniqFs = unique(GRID_SRATE_TARGET);
cache = struct();

for u = 1:numel(uniqFs)
    FsT = uniqFs(u);
    key = fskey(FsT);
    cache.(key) = cell(numel(files),1);

    if VERBOSE_STATUS
        fprintf('[CACHE] FsT=%.3f | %d arquivos...\n', FsT, numel(files));
    end

    for i = 1:numel(files)
        if VERBOSE_STATUS && mod(i,10)==0
            fprintf('  ... %d/%d\n', i, numel(files));
        end
        cache.(key){i} = preprocess_per_fs(raw{i}, FsT, VERBOSE_STATUS);
    end
end

%% ===================== ENUMERA CONFIGS (POUPANDO PARAMS QUANDO OFF) ===================== %%
cfg = enumerate_configs_gate_bp_harm_flow( ...
    GRID_INPUT_MODE, GRID_SRATE_TARGET, ...
    GRID_WI, GRID_WF, GRID_VOICES, GRID_CWT_WAVE, ...
    GRID_GATE_ON, GRID_HR_MIN_BPM, GRID_HR_MAX_BPM, GRID_DELTA_DB, GRID_SMOOTH_BINS, GRID_GATE_FLOOR, ...
    GRID_BP_ON, GRID_BP_ORDER, GRID_BP_WI, GRID_BP_WF, ...
    GRID_FLOW, ...
    GRID_BAND_W_HZ, GRID_DB_FLOOR_BAND, GRID_RIDGE_LAMBDA, GRID_F_JUMP_HZ, GRID_BURNIN_SEC, GRID_K_EVENT, GRID_PCHAIN_GAMMA, GRID_PCHAIN_MAX, ...
    GRID_HARM_ON, GRID_HARM_REL_DB, GRID_HARM_MIN_BPM, GRID_HARM_MAX_BPM, ...
    GRID_ECG_MIN_DIST_SEC, GRID_ECG_MIN_PROM, GRID_HR_MIN_BPM_ECG, GRID_HR_MAX_BPM_ECG, ...
    GRID_WIN_MEAN_SEC);

fprintf('Total configs: %d\n', numel(cfg));

%% ===================== PARALLEL SETUP (PROCESSES) ===================== %%
delete(gcp('nocreate'));
try
    parpool('Processes', MAX_WORKERS_SAFE);
catch
    parpool('local', MAX_WORKERS_SAFE);
end

rawConst   = parallel.pool.Constant(raw);
cacheConst = parallel.pool.Constant(cache);

if PAR_STATUS
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(ic) progress_print(ic, numel(cfg), PAR_PRINT_EVERY));
else
    dq = [];
end

%% ===================== RODA GRID (PARFOR) ===================== %%
tStart = tic;
rows = cell(numel(cfg),1);

parfor ic = 1:numel(cfg)
    C = cfg(ic);

    rmse_raw = []; mae_raw = []; corr_raw = [];
    rmse_sm  = []; mae_sm  = []; corr_sm  = [];
    n_ok = 0;

    n_ok_t     = zeros(nTypes,1);
    rmse_raw_t = cell(nTypes,1); mae_raw_t = cell(nTypes,1); corr_raw_t = cell(nTypes,1);
    rmse_sm_t  = cell(nTypes,1); mae_sm_t  = cell(nTypes,1); corr_sm_t  = cell(nTypes,1);

    rr = rawConst.Value;
    cc_struct = cacheConst.Value;

    key = fskey(C.SRATE_TARGET);
    if ~isfield(cc_struct, key)
        rows{ic} = config_row_types(C, ic, n_ok, rmse_raw, mae_raw, corr_raw, rmse_sm, mae_sm, corr_sm, ...
            typeTag, n_ok_t, rmse_raw_t, mae_raw_t, corr_raw_t, rmse_sm_t, mae_sm_t, corr_sm_t);
        if PAR_STATUS, send(dq, ic); end
        continue;
    end

    cc = cc_struct.(key);

    for fi = 1:numel(cc)
        try
            out = run_one_cached_gate_bp_harm(rr{fi}, cc{fi}, C);

            if out.has_ref
                n_ok = n_ok + 1;

                rmse_raw(end+1,1) = out.RMSE_raw; %#ok<PFOUS>
                mae_raw(end+1,1)  = out.MAE_raw;
                corr_raw(end+1,1) = out.CORR_raw;

                rmse_sm(end+1,1)  = out.RMSE_sm;
                mae_sm(end+1,1)   = out.MAE_sm;
                corr_sm(end+1,1)  = out.CORR_sm;

                ti = type_idx(fi);
                n_ok_t(ti) = n_ok_t(ti) + 1;

                rmse_raw_t{ti}(end+1,1) = out.RMSE_raw;
                mae_raw_t{ti}(end+1,1)  = out.MAE_raw;
                corr_raw_t{ti}(end+1,1) = out.CORR_raw;

                rmse_sm_t{ti}(end+1,1)  = out.RMSE_sm;
                mae_sm_t{ti}(end+1,1)   = out.MAE_sm;
                corr_sm_t{ti}(end+1,1)  = out.CORR_sm;
            end
        catch
        end
    end

    rows{ic} = config_row_types(C, ic, n_ok, rmse_raw, mae_raw, corr_raw, rmse_sm, mae_sm, corr_sm, ...
        typeTag, n_ok_t, rmse_raw_t, mae_raw_t, corr_raw_t, rmse_sm_t, mae_sm_t, corr_sm_t);

    if PAR_STATUS, send(dq, ic); end
end

Rall = vertcat(rows{:});
T = struct2table(Rall);

%% ===================== TABELAS FINAIS ===================== %%
T = sortrows(T, 'SCORE', 'ascend');

fprintf('\nConcluído em %.1fs\n', toc(tStart));
disp('=== TOP-20 GERAL (SCORE) ===');
disp(T(1:min(20,height(T)), :));

T_rmse = sortrows(T, 'RMSE_sm', 'ascend');
T_mae  = sortrows(T, 'MAE_sm',  'ascend');
T_corr = sortrows(T, 'CORR_sm', 'descend');

disp('=== TOP-10 por menor RMSE_sm ==='); disp(T_rmse(1:min(10,height(T_rmse)), :));
disp('=== TOP-10 por menor MAE_sm  ==='); disp(T_mae(1:min(10,height(T_mae)), :));
disp('=== TOP-10 por maior CORR_sm ==='); disp(T_corr(1:min(10,height(T_corr)), :));

for k = 1:nTypes
    colScore = "SCORE_" + typeTag(k);
    colNok   = "n_ok_"  + typeTag(k);

    if any(strcmp(T.Properties.VariableNames, colScore)) && any(strcmp(T.Properties.VariableNames, colNok))
        Tk = T(T.(colNok) > 0, :);
        Tk = sortrows(Tk, colScore, 'ascend');

        fprintf('\n=== TOP-10 (%s) por %s ===\n', typeTag(k), colScore);
        disp(Tk(1:min(10,height(Tk)), :));
    end
end

%% ===================== SALVA CSVs ===================== %%
out_all = fullfile(BASE_DIR, 'grid_GATE_BP_HARM_FLOW_parallel_ALL.csv');
writetable(T, out_all);
fprintf('\nSalvo (ALL): %s\n', out_all);

for k = 1:nTypes
    colScore = "SCORE_" + typeTag(k);
    colNok   = "n_ok_"  + typeTag(k);

    if any(strcmp(T.Properties.VariableNames, colScore))
        Tk = T(T.(colNok) > 0, :);
        Tk = sortrows(Tk, colScore, 'ascend');
        Tk = Tk(1:min(200,height(Tk)), :);

        out_k = fullfile(BASE_DIR, "grid_TOP_" + typeTag(k) + ".csv");
        writetable(Tk, out_k);
        fprintf('Salvo (TOP %s): %s\n', typeTag(k), out_k);
    end
end

%% ===================== FUNÇÕES ===================== %%

function progress_print(ic, n, every)
    if mod(ic, every)==0 || ic==1 || ic==n
        fprintf('[PAR] %d/%d (%.1f%%)\n', ic, n, 100*ic/n);
    end
end

function t = parse_type_num(nameNoExt)
    t = NaN;
    tok = regexp(nameNoExt, '_(\d+)(?:_|$)', 'tokens', 'once');
    if ~isempty(tok)
        t = str2double(tok{1});
    end
end

function S = load_one_mat(path)
    A = load(path);
    S = struct();
    S.path = path;

    S.has_radar = isfield(A,'radar_i') && isfield(A,'radar_q') && isfield(A,'fs_radar');
    if S.has_radar
        S.radar_i  = double(A.radar_i(:));
        S.radar_q  = double(A.radar_q(:));
        S.fs_radar = double(A.fs_radar);
    end

    S.has_ecg = isfield(A,'tfm_ecg1') && isfield(A,'fs_ecg');
    if S.has_ecg
        S.ecg1   = double(A.tfm_ecg1(:));
        S.fs_ecg = double(A.fs_ecg);
        if isfield(A,'tfm_ecg2'), S.ecg2 = double(A.tfm_ecg2(:)); else, S.ecg2 = []; end
    else
        S.ecg1 = []; S.ecg2 = []; S.fs_ecg = NaN;
    end
end

function key = fskey(Fs)
    key = sprintf('Fs_%g', Fs);
    key = strrep(key,'.','p');
end

function P = preprocess_per_fs(raw, FsT, VERBOSE_STATUS)
    P = struct();
    P.Fs = FsT;
    P.ok = true;

    if ~isfinite(FsT) || FsT <= 0
        P.ok = false; return;
    end

    % -------- radar -> fase -> resample --------
    if raw.has_radar && isfinite(raw.fs_radar) && raw.fs_radar > 0
        I   = force_finite_vector(raw.radar_i);
        Q   = force_finite_vector(raw.radar_q);
        fs0 = raw.fs_radar;

        phase = unwrap(atan2(Q, I));
        x0    = wrap_phase(phase);

        t0 = (0:numel(x0)-1)'/fs0;

        [p,q] = rat(FsT/fs0, 1e-12);
        x  = resample(x0, p, q);
        Fs = fs0 * p / q;

        x = force_finite_vector(x);
        t = (0:numel(x)-1)'/Fs;
        if ~isempty(t0) && ~isempty(t) && t(end) > t0(end)
            t(end) = t0(end);
        end

        P.radar_ok = true;
        P.radar_t  = t;
        P.radar_x  = x(:);
        P.radar_Fs = Fs;
    else
        P.radar_ok = false;
        P.radar_t  = [];
        P.radar_x  = [];
        P.radar_Fs = NaN;
    end

    % -------- ECG pré-processado p/ referência --------
    P.ecg_ok  = false;
    P.ecg_fs  = NaN;
    P.ecg_t   = [];
    P.ecg_f   = [];
    P.ecg_mwi = [];

    if raw.has_ecg && isfinite(raw.fs_ecg) && raw.fs_ecg > 0 && ~isempty(raw.ecg1)
        fs  = raw.fs_ecg;
        ecg = force_finite_vector(raw.ecg1);

        try
            w_baseline = max(3, round(0.20 * fs));
            ecg0 = ecg - movmedian(ecg, w_baseline, 'Endpoints','shrink');
            ecg0 = force_finite_vector(ecg0);

            f1 = 5;  f2 = 25;
            if f2 >= fs/2, f2 = 0.9*(fs/2); end
            if f1 <= 0 || f1 >= f2
                if VERBOSE_STATUS
                    fprintf('[ECG PRE] cortes inválidos (%s)\n', raw.path);
                end
            else
                [b,a] = butter(2, [f1 f2]/(fs/2), 'bandpass');
                ecg_f = filtfilt(b,a, ecg0);
                ecg_f = force_finite_vector(ecg_f);

                ecg_f = ecg_f / (mad(ecg_f,1) + eps);

                d  = [0; diff(ecg_f)];
                d  = d / (mad(d,1) + eps);
                sq = d.^2;

                w_int = max(3, round(0.12 * fs));
                mwi = movmean(sq, w_int, 'Endpoints','shrink');
                mwi = force_finite_vector(mwi);

                t_ecg = (0:numel(ecg_f)-1)'/fs;

                P.ecg_ok  = (numel(ecg_f) >= 10) && all(isfinite(ecg_f)) && all(isfinite(mwi));
                P.ecg_fs  = fs;
                P.ecg_t   = t_ecg;
                P.ecg_f   = ecg_f;
                P.ecg_mwi = mwi;
            end
        catch
            if VERBOSE_STATUS
                fprintf('[ECG PRE FAIL] %s\n', raw.path);
            end
            P.ecg_ok = false;
            P.ecg_fs = fs;
        end
    end
end

function out = run_one_cached_gate_bp_harm(raw, pre, C)
    out = struct('has_ref',false,'RMSE_raw',NaN,'MAE_raw',NaN,'CORR_raw',NaN,'RMSE_sm',NaN,'MAE_sm',NaN,'CORR_sm',NaN);

    % ------- seleciona entrada -------
    switch lower(string(C.INPUT_MODE))
        case "radar"
            if ~pre.radar_ok, return; end
            x_puro = pre.radar_x;
            t      = pre.radar_t;
            Fs     = pre.radar_Fs;
        otherwise
            return;
    end

    x_puro = force_finite_vector(x_puro);
    if ~isfinite(Fs) || Fs <= 0, return; end

    % ------- validações básicas -------
    if C.WF >= Fs/2, return; end
    if C.WI <= 0 || C.WI >= C.WF, return; end

    % ------- aplica pipeline (comparação de ordem) -------
    x_fin = x_puro;

    % Se GATE_OFF, força FLOW='NA' (não compara) — mas mantém campo no output
    if C.GATE_ON ~= 1
        % Gate OFF: faz só BP (se existir)
        if C.BP_ON == 1
            x_fin = apply_bp_sos(x_fin, Fs, C.BP_ORDER, C.BP_WI, C.BP_WF);
        end
    else
        % Gate ON: depende do FLOW
        if strcmpi(C.FLOW,'BP_THEN_GATE')
            if C.BP_ON == 1
                x_fin = apply_bp_sos(x_fin, Fs, C.BP_ORDER, C.BP_WI, C.BP_WF);
            end
            x_fin = gate_by_power_hr_full(x_fin, Fs, C.HR_MIN_BPM, C.HR_MAX_BPM, C.DELTA_DB, C.SMOOTH_BINS, C.GATE_FLOOR);
        else
            % default: GATE_THEN_BP
            x_fin = gate_by_power_hr_full(x_fin, Fs, C.HR_MIN_BPM, C.HR_MAX_BPM, C.DELTA_DB, C.SMOOTH_BINS, C.GATE_FLOOR);
            if C.BP_ON == 1
                x_fin = apply_bp_sos(x_fin, Fs, C.BP_ORDER, C.BP_WI, C.BP_WF);
            end
        end
    end

    x_fin = force_finite_vector(x_fin);

    % ------- CWT -------
    [wt, fbin] = cwt(x_fin, char(C.CWT_WAVE), Fs, ...
        'FrequencyLimits',[C.WI C.WF], ...
        'VoicesPerOctave', C.VOICES);

    fbin = fbin(:);
    if numel(fbin)>1 && fbin(2)<fbin(1)
        fbin = flipud(fbin);
        wt   = flipud(wt);
    end

    P = single(abs(wt).^2);
    clear wt;

    if any(~isfinite(P(:)))
        P(~isfinite(P)) = 0;
    end

    % ------- Ridge (área + punish + evento + harm) -------
    freq_hz = ridge_area_punish_event_harm(P, fbin, Fs, ...
        C.BAND_W_HZ, C.DB_FLOOR_BAND, ...
        C.RIDGE_LAMBDA, C.F_JUMP_HZ, C.BURNIN_SEC, C.K_EVENT, C.PCHAIN_GAMMA, C.PCHAIN_MAX, ...
        C.HARM_ON, C.HARM_REL_DB, C.HARM_MIN_BPM, C.HARM_MAX_BPM);

    est_bpm = 60*freq_hz(:);

    % ------- ECG ref -------
    hr_ref = hr_ref_on_t(pre, t, C.ECG_MIN_DIST_SEC, C.ECG_MIN_PROM, C.HR_MIN_BPM_ECG, C.HR_MAX_BPM_ECG);
    if isempty(hr_ref) || all(~isfinite(hr_ref)), return; end

    okR = isfinite(est_bpm) & isfinite(hr_ref);
    if nnz(okR) < 10, return; end

    eR = est_bpm(okR) - hr_ref(okR);
    out.has_ref  = true;
    out.RMSE_raw = sqrt(mean(eR.^2,'omitnan'));
    out.MAE_raw  = mean(abs(eR),'omitnan');
    out.CORR_raw = corr(est_bpm(okR), hr_ref(okR), 'Rows','complete');

    % ------- movmean -------
    winS   = max(1, round(C.WIN_MEAN_SEC * Fs));
    est_sm = movmean(est_bpm, winS, 'Endpoints','shrink');
    ref_sm = movmean(hr_ref,  winS, 'Endpoints','shrink');

    okS = isfinite(est_sm) & isfinite(ref_sm);
    if nnz(okS) < 10, return; end

    eS = est_sm(okS) - ref_sm(okS);
    out.RMSE_sm = sqrt(mean(eS.^2,'omitnan'));
    out.MAE_sm  = mean(abs(eS),'omitnan');
    out.CORR_sm = corr(est_sm(okS), ref_sm(okS), 'Rows','complete');
end

function x = apply_bp_sos(x, Fs, ord, wi, wf)
    x = force_finite_vector(x);

    if ~isfinite(Fs) || Fs<=0, return; end
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

function hr = hr_ref_on_t(pre, t, minDistSec, minProm, hrMin, hrMax)
    hr = [];
    if ~pre.ecg_ok, return; end

    fs    = pre.ecg_fs;
    mwi   = pre.ecg_mwi(:);
    ecg_f = pre.ecg_f(:);
    t_ecg = pre.ecg_t(:);

    if numel(mwi) < 10 || numel(ecg_f) ~= numel(mwi), return; end

    thr = median(mwi) + 3.0*mad(mwi,1);

    minDistSec2 = max(minDistSec, 60/max(hrMax,eps));
    minDistSmp  = max(1, round(minDistSec2 * fs));

    prom_adapt = max(minProm, 1.5*mad(mwi,1));

    try
        [~, locs_env] = findpeaks(mwi, ...
            'MinPeakDistance',  minDistSmp, ...
            'MinPeakHeight',    thr, ...
            'MinPeakProminence', prom_adapt);
    catch
        return;
    end

    if numel(locs_env) < 3, return; end

    w_ref = max(1, round(0.08 * fs));
    n = numel(ecg_f);
    locs_R = zeros(size(locs_env));

    for k = 1:numel(locs_env)
        c = locs_env(k);
        i1 = max(1, c - w_ref);
        i2 = min(n, c + w_ref);
        [~, im] = max(abs(ecg_f(i1:i2)));
        locs_R(k) = i1 + im - 1;
    end

    locs_R = unique(locs_R,'stable');
    if numel(locs_R) < 3, return; end

    t_R = t_ecg(locs_R);
    RR = diff(t_R);
    if isempty(RR), return; end

    hr_inst = 60 ./ RR;
    t_hr    = t_R(1:end-1) + RR/2;

    ok = isfinite(hr_inst) & isfinite(t_hr) & (hr_inst >= hrMin) & (hr_inst <= hrMax);
    hr_inst = hr_inst(ok);
    t_hr    = t_hr(ok);

    if numel(hr_inst) < 3, return; end

    hr = interp1(t_hr, hr_inst, t, 'linear', 'extrap');
end

function R = config_row_types(C, cfg_id, n_ok, rmse_raw, mae_raw, corr_raw, rmse_sm, mae_sm, corr_sm, ...
                             typeTag, n_ok_t, rmse_raw_t, mae_raw_t, corr_raw_t, rmse_sm_t, mae_sm_t, corr_sm_t)
    R = C;
    R.cfg_id = cfg_id;

    R.n_ok     = n_ok;
    R.RMSE_raw = mean(rmse_raw,'omitnan');
    R.MAE_raw  = mean(mae_raw,'omitnan');
    R.CORR_raw = mean(corr_raw,'omitnan');

    R.RMSE_sm  = mean(rmse_sm,'omitnan');
    R.MAE_sm   = mean(mae_sm,'omitnan');
    R.CORR_sm  = mean(corr_sm,'omitnan');

    R.SCORE = R.RMSE_sm - 5*R.CORR_sm;

    for k = 1:numel(typeTag)
        tag = char(typeTag(k));

        R.(['n_ok_' tag])     = n_ok_t(k);

        R.(['RMSE_raw_' tag]) = mean(rmse_raw_t{k},'omitnan');
        R.(['MAE_raw_'  tag]) = mean(mae_raw_t{k}, 'omitnan');
        R.(['CORR_raw_' tag]) = mean(corr_raw_t{k},'omitnan');

        R.(['RMSE_sm_'  tag]) = mean(rmse_sm_t{k}, 'omitnan');
        R.(['MAE_sm_'   tag]) = mean(mae_sm_t{k},  'omitnan');
        R.(['CORR_sm_'  tag]) = mean(corr_sm_t{k}, 'omitnan');

        R.(['SCORE_' tag])    = R.(['RMSE_sm_' tag]) - 5*R.(['CORR_sm_' tag]);
    end
end

function cfg = enumerate_configs_gate_bp_harm_flow( ...
    GRID_INPUT_MODE, GRID_SRATE_TARGET, ...
    GRID_WI, GRID_WF, GRID_VOICES, GRID_CWT_WAVE, ...
    GRID_GATE_ON, GRID_HR_MIN_BPM, GRID_HR_MAX_BPM, GRID_DELTA_DB, GRID_SMOOTH_BINS, GRID_GATE_FLOOR, ...
    GRID_BP_ON, GRID_BP_ORDER, GRID_BP_WI, GRID_BP_WF, ...
    GRID_FLOW, ...
    GRID_BAND_W_HZ, GRID_DB_FLOOR_BAND, GRID_RIDGE_LAMBDA, GRID_F_JUMP_HZ, GRID_BURNIN_SEC, GRID_K_EVENT, GRID_PCHAIN_GAMMA, GRID_PCHAIN_MAX, ...
    GRID_HARM_ON, GRID_HARM_REL_DB, GRID_HARM_MIN_BPM, GRID_HARM_MAX_BPM, ...
    GRID_ECG_MIN_DIST_SEC, GRID_ECG_MIN_PROM, GRID_HR_MIN_BPM_ECG, GRID_HR_MAX_BPM_ECG, ...
    GRID_WIN_MEAN_SEC)

    cfg = struct([]);
    ic = 0;

    % garante FLOW como cellstr
    if isstring(GRID_FLOW), GRID_FLOW = cellstr(GRID_FLOW); end
    if isempty(GRID_FLOW),  GRID_FLOW = {'GATE_THEN_BP'};   end

    for INPUT_MODE = GRID_INPUT_MODE
    for SRATE_TARGET = GRID_SRATE_TARGET
    for WI = GRID_WI
    for WF = GRID_WF
    for VOICES = GRID_VOICES
    for CWT_WAVE = GRID_CWT_WAVE

    for BP_ON = GRID_BP_ON
    for BP_ORDER = GRID_BP_ORDER
    for BP_WI = GRID_BP_WI
    for BP_WF = GRID_BP_WF

    for GATE_ON = GRID_GATE_ON

        % se gate desligado, poupa params do gate e a comparação de FLOW
        if GATE_ON == 1
            DELTA_LIST  = GRID_DELTA_DB;      % numérico
            SMOOTH_LIST = GRID_SMOOTH_BINS;   % numérico
            FLOOR_LIST  = GRID_GATE_FLOOR;    % numérico
            FLOW_LIST   = GRID_FLOW;          % cellstr
        else
            DELTA_LIST  = NaN;
            SMOOTH_LIST = NaN;
            FLOOR_LIST  = NaN;
            FLOW_LIST   = {'NA'};
        end

        for HR_MIN_BPM = GRID_HR_MIN_BPM
        for HR_MAX_BPM = GRID_HR_MAX_BPM
        for DELTA_DB = DELTA_LIST
        for SMOOTH_BINS = SMOOTH_LIST
        for GATE_FLOOR = FLOOR_LIST
        for FLOW = FLOW_LIST

            % se BP desligado (se você adicionar BP_ON=0), poupa params do BP
            if BP_ON ~= 1
                BP_ORDER2 = NaN; BP_WI2 = NaN; BP_WF2 = NaN;
            else
                BP_ORDER2 = BP_ORDER; BP_WI2 = BP_WI; BP_WF2 = BP_WF;
            end

            for BAND_W_HZ = GRID_BAND_W_HZ
            for DB_FLOOR_BAND = GRID_DB_FLOOR_BAND
            for RIDGE_LAMBDA = GRID_RIDGE_LAMBDA
            for F_JUMP_HZ = GRID_F_JUMP_HZ
            for BURNIN_SEC = GRID_BURNIN_SEC
            for K_EVENT = GRID_K_EVENT
            for PCHAIN_GAMMA = GRID_PCHAIN_GAMMA
            for PCHAIN_MAX = GRID_PCHAIN_MAX

            for HARM_ON = GRID_HARM_ON

                % se harm desligado, poupa params do harm
                if HARM_ON == 1
                    HREL_LIST = GRID_HARM_REL_DB;    % numérico
                    HMIN_LIST = GRID_HARM_MIN_BPM;   % numérico
                    HMAX_LIST = GRID_HARM_MAX_BPM;   % numérico
                else
                    HREL_LIST = NaN;
                    HMIN_LIST = NaN;
                    HMAX_LIST = NaN;
                end

                for HARM_REL_DB = HREL_LIST
                for HARM_MIN_BPM = HMIN_LIST
                for HARM_MAX_BPM = HMAX_LIST

                    for ECG_MIN_DIST_SEC = GRID_ECG_MIN_DIST_SEC
                    for ECG_MIN_PROM = GRID_ECG_MIN_PROM
                    for HR_MIN_BPM_ECG = GRID_HR_MIN_BPM_ECG
                    for HR_MAX_BPM_ECG = GRID_HR_MAX_BPM_ECG
                    for WIN_MEAN_SEC = GRID_WIN_MEAN_SEC

                        % validações
                        if WI <= 0 || WI >= WF, continue; end
                        if (HR_MIN_BPM >= HR_MAX_BPM), continue; end
                        if (HR_MIN_BPM_ECG >= HR_MAX_BPM_ECG), continue; end
                        if BP_ON == 1 && (BP_WI2 <= 0 || BP_WI2 >= BP_WF2), continue; end

                        % compara FLOW só quando gate=1 e bp=1
                        if (GATE_ON==1) && (BP_ON~=1)
                            FLOW = {'NA'};
                        end

                        ic = ic + 1;

                        % sempre preencher TODOS os campos (evita "dissimilar structures")
                        cfg(ic).INPUT_MODE    = char(INPUT_MODE);
                        cfg(ic).SRATE_TARGET  = SRATE_TARGET;

                        cfg(ic).WI       = WI;
                        cfg(ic).WF       = WF;
                        cfg(ic).VOICES   = VOICES;
                        cfg(ic).CWT_WAVE = char(CWT_WAVE);

                        cfg(ic).GATE_ON      = GATE_ON;
                        cfg(ic).HR_MIN_BPM   = HR_MIN_BPM;
                        cfg(ic).HR_MAX_BPM   = HR_MAX_BPM;
                        cfg(ic).DELTA_DB     = DELTA_DB;
                        cfg(ic).SMOOTH_BINS  = SMOOTH_BINS;
                        cfg(ic).GATE_FLOOR   = GATE_FLOOR;

                        cfg(ic).BP_ON    = BP_ON;
                        cfg(ic).BP_ORDER = BP_ORDER2;
                        cfg(ic).BP_WI    = BP_WI2;
                        cfg(ic).BP_WF    = BP_WF2;

                        cfg(ic).FLOW = char(FLOW);  % FLOW é cellstr => char ok

                        cfg(ic).BAND_W_HZ      = BAND_W_HZ;
                        cfg(ic).DB_FLOOR_BAND  = DB_FLOOR_BAND;
                        cfg(ic).RIDGE_LAMBDA   = RIDGE_LAMBDA;
                        cfg(ic).F_JUMP_HZ      = F_JUMP_HZ;
                        cfg(ic).BURNIN_SEC     = BURNIN_SEC;
                        cfg(ic).K_EVENT        = K_EVENT;
                        cfg(ic).PCHAIN_GAMMA   = PCHAIN_GAMMA;
                        cfg(ic).PCHAIN_MAX     = PCHAIN_MAX;

                        cfg(ic).HARM_ON       = HARM_ON;
                        cfg(ic).HARM_REL_DB   = HARM_REL_DB;
                        cfg(ic).HARM_MIN_BPM  = HARM_MIN_BPM;
                        cfg(ic).HARM_MAX_BPM  = HARM_MAX_BPM;

                        cfg(ic).ECG_MIN_DIST_SEC = ECG_MIN_DIST_SEC;
                        cfg(ic).ECG_MIN_PROM     = ECG_MIN_PROM;
                        cfg(ic).HR_MIN_BPM_ECG   = HR_MIN_BPM_ECG;
                        cfg(ic).HR_MAX_BPM_ECG   = HR_MAX_BPM_ECG;

                        cfg(ic).WIN_MEAN_SEC = WIN_MEAN_SEC;

                    end; end; end; end; end
                end; end; end
            end; end; end; end; end; end; end; end
        end; end; end; end; end; end
        end; end; end; end; end; end; end; end; end; end; end; end;
end


%% ===== GATE (robusto: broad + spike) =====
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

%% ===== RIDGE (área + punish + evento + harm f/2) =====
function freq_hz = ridge_area_punish_event_harm(P, f_bin, srate, BAND_W_HZ, DB_FLOOR_BAND, RIDGE_LAMBDA, F_JUMP_HZ, BURNIN_SEC, K_EVENT, GAMMA, KMAX, ...
                                                HARM_ON, HARM_REL_DB, HARM_MIN_BPM, HARM_MAX_BPM)
    f_bin = double(f_bin(:));
    nb = numel(f_bin);
    nt = size(P,2);

    df = gradient(f_bin);
    df(df<=0) = eps;

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

    % mapa "metade da freq" para um bin inicial viável
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

    % buffers (evita realocar no loop)
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

        % pico local por faixa (loop curto; faixa pequena)
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

        % ---- HARM f/2 ----
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
