% plot_csv_data.m
% Comparativos HR e RR
% - Bases HR/RR: CSV "dt_seconds,value" (ignora '#', separador ',')
% - Overlay HR (Polar): CSV com ';' e cabeçalho "Phone timestamp;HR [bpm];..."
% - Overlay HR/RR (OmniBuds): CSV com metadados; gravação inicia na linha 17;
%   usar SEMPRE as colunas 3 (valor) e 4 (timestamp, possivelmente 'quoted').
% - t0 do segundo dataset (Omni) = primeiro timestamp que aparecer (entre HR e RR).
% Saídas: HR_compare.png e RR_compare.png

clc; clear; close all;

% ========================
% === ARQUIVOS / FLAGS ===
% ========================
% Bases (dt_seconds,value)
base_hr_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes.csv';      % base HR
base_rr_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes_rr.csv';   % base RR

% Overlay HR (Polar, separado por ';')
polar_hr_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLAR2.csv';  % overlay HR (Polar)

% Overlay no segundo formato (OmniBuds)
omni_hr_path  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\omni\Heart1.csv';  % overlay HR (Omni)
omni_rr_path  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\omni\Breath1.csv'; % overlay RR (Omni)

% Escolha da fonte de HR overlay
USE_OMNI_FOR_HR = true;   % true => usa OmniBuds para HR, false => usa Polar

% ====== Config de erro ======
RESCALE_OVERLAY_BEFORE_ERROR = true;   % reescalar overlay para a faixa do base antes do erro
RESCALE_METHOD = 'range';              % 'range' ou 'zscore'

% ====== Janelas de média móvel ======
MOV_WIN_BASE = 15;   % amostras (para base)
MOV_WIN_OVLY = 15;   % amostras (para overlays)

% =========================
% === LEITURA DAS BASES ===
% =========================
fprintf('Lendo base HR: %s\n', base_hr_path);
[t_hr_base, v_hr_base, v_hr_base_mov] = read_base_dt_val(base_hr_path, MOV_WIN_BASE);

fprintf('Lendo base RR: %s\n', base_rr_path);
[t_rr_base, v_rr_base, v_rr_base_mov] = read_base_dt_val(base_rr_path, MOV_WIN_BASE);

% ==================================================
% === DETECTAR t0 ABSOLUTO COMPARTILHADO (OMNI)  ===
% ==================================================
% Pré-ler apenas o t0 absoluto dos arquivos Omni (se existirem),
% para definir um t0 comum entre HR e RR do segundo dataset.
omni_t0_abs = NaT;

if exist(omni_hr_path,'file') == 2
    hr_t0_abs = peek_omni_first_dt(omni_hr_path);
    if ~isnat(hr_t0_abs), omni_t0_abs = hr_t0_abs; end
end
if exist(omni_rr_path,'file') == 2
    rr_t0_abs = peek_omni_first_dt(omni_rr_path);
    if ~isnat(rr_t0_abs)
        if isnat(omni_t0_abs)
            omni_t0_abs = rr_t0_abs;
        else
            omni_t0_abs = min(omni_t0_abs, rr_t0_abs);
        end
    end
end

% ===================================
% === 1) HR: base mov vs overlays ===
% ===================================
hr_t = []; hr_v = []; hr_v_mov = [];

if USE_OMNI_FOR_HR && exist(omni_hr_path,'file') == 2
    fprintf('Lendo HR (OmniBuds): %s\n', omni_hr_path);
    [hr_t_rel_from_omni_t0, hr_v, hr_v_mov] = read_omni_value_timestamp(omni_hr_path, MOV_WIN_OVLY, omni_t0_abs);
    if ~isempty(hr_t_rel_from_omni_t0)
        % Alinhar o overlay ao início da base HR mantendo a escala relativa do Omni (t0 comum)
        hr_t = hr_t_rel_from_omni_t0 + t_hr_base(1);
    end
else
    fprintf('Lendo HR (Polar): %s\n', polar_hr_path);
    [hr_t, hr_v, hr_v_mov] = read_polar_hr(polar_hr_path, t_hr_base(1), MOV_WIN_OVLY);
end

% === Figura HR (sem traço bruto da base) ===
figHR = figure('Name','HR: Base (movmean) vs Overlay','Color','w');
plot(t_hr_base, v_hr_base_mov, 'r', 'LineWidth', 1.3, ...
     'DisplayName', sprintf('HR base mov(%d)', MOV_WIN_BASE)); hold on;

if ~isempty(hr_t)
    plot(hr_t, hr_v, 'g', 'LineWidth', 1.0, 'DisplayName', 'HR overlay (bruto)');
    if ~isempty(hr_v_mov)
        plot(hr_t, hr_v_mov, '--k', 'LineWidth', 1.1, ...
             'DisplayName', sprintf('HR overlay mov(%d)', MOV_WIN_OVLY));
    end
end

grid on; xlabel('Tempo (s)'); ylabel('HR [bpm]');
title('HR — Base (movmean) vs Overlay');
[xmin_hr, xmax_hr, ymin_hr, ymax_hr] = xy_limits_union( ...
    t_hr_base, v_hr_base_mov, hr_t, hr_v, hr_v_mov);
xlim([xmin_hr xmax_hr]); ylim([ymin_hr ymax_hr]);
legend('Location','best');

% Métricas HR (base mov vs overlay bruto)
if ~isempty(hr_t)
    [rmse_hr, mae_hr, mape_hr] = error_metrics(t_hr_base, v_hr_base_mov, hr_t, hr_v, ...
        RESCALE_OVERLAY_BEFORE_ERROR, RESCALE_METHOD);
    fprintf('\n=== ERRO HR (base mov vs overlay%s) ===\n', ...
        ternary(RESCALE_OVERLAY_BEFORE_ERROR, ' [reescalado]', ' [bruto]'));
    fprintf('RMSE: %.4f\n', rmse_hr);
    fprintf('MAE : %.4f\n', mae_hr);
    fprintf('MAPE: %.2f %%\n', mape_hr);
else
    fprintf('\nSem overlay de HR para cálculo de erro.\n');
end

exportgraphics(figHR, 'HR_compare.png','Resolution',150);

% ==================================
% === 2) RR: base mov vs overlays ===
% ==================================
rr_t = []; rr_v = []; rr_v_mov = [];

if exist(omni_rr_path,'file') == 2
    fprintf('Lendo RR (OmniBuds): %s\n', omni_rr_path);
    [rr_t_rel_from_omni_t0, rr_v, rr_v_mov] = read_omni_value_timestamp(omni_rr_path, MOV_WIN_OVLY, omni_t0_abs);
    if ~isempty(rr_t_rel_from_omni_t0)
        % Alinhar o overlay ao início da base RR mantendo a escala relativa do Omni (t0 comum)
        rr_t = rr_t_rel_from_omni_t0 + t_rr_base(1);
    end
else
    fprintf('Overlay RR (OmniBuds) não encontrado: %s\n', omni_rr_path);
end

% === Figura RR (sem traço bruto da base) ===
figRR = figure('Name','RR: Base (movmean) vs Overlay','Color','w');
plot(t_rr_base, v_rr_base_mov, 'r', 'LineWidth', 1.3, ...
     'DisplayName', sprintf('RR base mov(%d)', MOV_WIN_BASE)); hold on;


    if ~isempty(rr_v_mov)
        plot(rr_t, rr_v_mov, '--k', 'LineWidth', 1.1, ...
             'DisplayName', sprintf('RR overlay mov(%d)', MOV_WIN_OVLY));
    end

grid on; xlabel('Tempo (s)'); ylabel('RR [rpm]');
title('RR — Base (movmean) vs Overlay');
[xmin_rr, xmax_rr, ymin_rr, ymax_rr] = xy_limits_union( ...
    t_rr_base, v_rr_base_mov, rr_t, rr_v, rr_v_mov);
xlim([xmin_rr xmax_rr]); ylim([ymin_rr ymax_rr]);
legend('Location','best');

% Métricas RR (base mov vs overlay bruto)
if ~isempty(rr_t)
    [rmse_rr, mae_rr, mape_rr] = error_metrics(t_rr_base, v_rr_base_mov, rr_t, rr_v, ...
        RESCALE_OVERLAY_BEFORE_ERROR, RESCALE_METHOD);
    fprintf('\n=== ERRO RR (base mov vs overlay%s) ===\n', ...
        ternary(RESCALE_OVERLAY_BEFORE_ERROR, ' [reescalado]', ' [bruto]'));
    fprintf('RMSE: %.4f\n', rmse_rr);
    fprintf('MAE : %.4f\n', mae_rr);
    fprintf('MAPE: %.2f %%\n', mape_rr);
else
    fprintf('\nSem overlay de RR para cálculo de erro.\n');
end

exportgraphics(figRR, 'RR_compare.png','Resolution',150);

% ==================
% === FUNÇÕES   ===
% ==================

function [t, v, v_mov] = read_base_dt_val(p, mov_win)
    if exist(p,'file') ~= 2
        error('Arquivo base não encontrado: %s', p);
    end
    data = readmatrix(p, 'CommentStyle','#', 'Delimiter', ',');
    if isempty(data) || size(data,2) < 2
        error('Base inválida (precisa "dt_seconds,value"): %s', p);
    end
    t = data(:,1);
    v = data(:,2);
    v_mov = movmean(v, mov_win, 'omitnan');
end

function [th_s, vh, vh_mov] = read_polar_hr(p, t0_align, mov_win)
    % Lê POLAR (;) com colunas "Phone timestamp" e "HR [bpm]"
    th_s = []; vh = []; vh_mov = [];
    if exist(p,'file') ~= 2, return; end

    % detectImportOptions com fallback (compatibilidade)
    try
        opts = detectImportOptions(p, ...
            'Delimiter',';', ...
            'ReadVariableNames',true, ...
            'TextType','string', ...
            'ConsecutiveDelimitersRule','join', ...
            'LeadingDelimitersRule','ignore');
        T = readtable(p, opts);
    catch
        % Fallback genérico
        T = readtable(p, ...
            'Delimiter',';', ...
            'ReadVariableNames',true, ...
            'TextType','string');
    end

    norm = lower(regexprep(T.Properties.VariableNames, '[^a-zA-Z0-9]', ''));
    idx_ts = find(strcmp(norm, 'phonetimestamp'), 1);
    idx_hr = find(strcmp(norm, 'hrbpm'), 1);
    if isempty(idx_ts), idx_ts = 1; end
    if isempty(idx_hr), idx_hr = min(2, width(T)); end

    ts_raw = string(T{:, idx_ts});
    hr_raw = string(T{:, idx_hr});
    vh = str2double(strrep(hr_raw, ',', '.'));

    th = NaT(size(ts_raw));
    ok = false;
    fmts = {'yyyy-MM-dd''T''HH:mm:ss.SSS','yyyy-MM-dd''T''HH:mm:ss', ...
            'dd/MM/yyyy HH:mm:ss.SSS','dd/MM/yyyy HH:mm:ss'};
    for k = 1:numel(fmts)
        try
            th = datetime(ts_raw, 'InputFormat', fmts{k}, 'TimeZone','local');
            ok = true; break;
        catch, end
    end
    if ~ok
        try, th = datetime(ts_raw, 'TimeZone','local'); catch, th = NaT(size(ts_raw)); end
    end

    mask = ~isnat(th) & ~isnan(vh);
    th = th(mask); vh = vh(mask);
    if isempty(th), return; end

    % segundos relativos e alinhamento em t0_align
    th_s = seconds(th - th(1));
    th_s = th_s + t0_align;
    vh_mov = movmean(vh, mov_win, 'omitnan');
end

function first_abs_dt = peek_omni_first_dt(p)
    % Lê APENAS o primeiro timestamp absoluto do arquivo Omni (linha 17, col 4)
    first_abs_dt = NaT;
    if exist(p,'file') ~= 2, return; end

    % Tentativa 1: detectImportOptions + DataLines (se disponível)
    try
        opts = detectImportOptions(p, 'Delimiter', ',', ...
            'ConsecutiveDelimitersRule','join', 'TextType','string', ...
            'ReadVariableNames', false);
        % Se a versão suportar, definimos DataLines
        try
            opts.DataLines = [17 17];         % pode não existir em versões antigas
            T1 = readtable(p, opts);
        catch
            % Tentativa 2: readtable com HeaderLines (compat antigo)
            T1 = readtable(p, 'Delimiter', ',', 'ReadVariableNames', false, ...
                'HeaderLines', 16, 'TextType','string', ...
                'ConsecutiveDelimitersRule','join', 'LeadingDelimitersRule','ignore');
            if height(T1) > 1
                T1 = T1(1,:); % só a 1ª linha após o cabeçalho
            end
        end
    catch
        % Tentativa 3: textscan
        T1 = omni_textscan_single_row(p);
    end

    if isempty(T1) || width(T1) < 4 || height(T1) == 0, return; end

    ts = string(T1{1,4});
    ts = erase(ts, '''');
    tsn = str2double(ts);
    if isnan(tsn), return; end

    % Segundos vs milissegundos
    if tsn > 1e11, secs = tsn/1000; else, secs = tsn; end
    try
        first_abs_dt = datetime(secs, 'ConvertFrom','posixtime', 'TimeZone','UTC');
    catch
        first_abs_dt = NaT;
    end
end

function [t_rel_from_global_t0, val, val_mov] = read_omni_value_timestamp(p, mov_win, omni_t0_abs)
    % Lê OmniBuds a partir da LINHA 17, usando colunas 3 (valor) e 4 (timestamp)
    % Compatível com versões antigas (sem DataLines).
    t_rel_from_global_t0 = []; val = []; val_mov = [];
    if exist(p,'file') ~= 2, return; end

    % Tentativa 1: detectImportOptions + DataLines
    used = false;
    try
        opts = detectImportOptions(p, 'Delimiter', ',', ...
            'ConsecutiveDelimitersRule','join', 'TextType','string', ...
            'ReadVariableNames', false);
        try
            opts.DataLines = [17 Inf];       % pode não existir em versões antigas
            T = readtable(p, opts);
            used = true;
        catch
            used = false;
        end
    catch
        used = false;
    end

    % Tentativa 2: readtable com HeaderLines (compat antigo)
    if ~used
        try
            T = readtable(p, ...
                'Delimiter', ',', ...
                'ReadVariableNames', false, ...
                'HeaderLines', 16, ...
                'TextType', 'string', ...
                'ConsecutiveDelimitersRule','join', ...
                'LeadingDelimitersRule','ignore');
            used = true;
        catch
            used = false;
        end
    end

    % Tentativa 3: textscan (fallback final)
    if ~used
        T = omni_textscan_all_rows(p);
        if isempty(T), return; end
    end

    if width(T) < 4 || height(T) == 0
        warning('OmniBuds: tabela vazia ou com colunas insuficientes: %s', p);
        return;
    end

    sval = string(T{:,3}); sval = erase(sval, ''''); sval = strrep(sval, ',', '.');
    sts  = string(T{:,4}); sts  = erase(sts,  '''');

    val = str2double(sval);
    tsn = str2double(sts);
    m = ~isnan(val) & ~isnan(tsn);
    val = val(m); tsn = tsn(m);
    if isempty(tsn), return; end

    % Detecta ms vs s
    secs = tsn;
    if any(tsn > 1e11), secs = tsn/1000; end

    % Datetimes absolutos
    try
        dt_abs = datetime(secs, 'ConvertFrom','posixtime', 'TimeZone','UTC');
    catch
        dt_abs = datetime(1970,1,1,'TimeZone','UTC') + seconds(secs);
    end

    % t0 absoluto comum
    if isnat(omni_t0_abs), omni_t0_abs = dt_abs(1); end

    % Tempo relativo ao t0 ABSOLUTO COMUM
    t_rel_from_global_t0 = seconds(dt_abs - omni_t0_abs);
    val_mov = movmean(val, mov_win, 'omitnan');
end

function T1 = omni_textscan_single_row(p)
    % Fallback: ler somente uma linha (a 1ª após 16 linhas de cabeçalho)
    T1 = [];
    fid = fopen(p,'r');
    if fid < 0, return; end
    cleaner = onCleanup(@() fclose(fid));
    % Pula 16 linhas
    for k=1:16
        if feof(fid), return; end
        fgetl(fid);
    end
    if feof(fid), return; end
    C = textscan(fid, '%s%s%s%s', 1, 'Delimiter', ',', ...
        'ReturnOnError', false, 'MultipleDelimsAsOne', true, ...
        'Whitespace','', 'EndOfLine','\n');
    if isempty(C) || numel(C) < 4, return; end
    T1 = table(string(C{1}), string(C{2}), string(C{3}), string(C{4}));
end

function T = omni_textscan_all_rows(p)
    % Fallback: ler todas as linhas a partir da 17ª
    T = [];
    fid = fopen(p,'r');
    if fid < 0, return; end
    cleaner = onCleanup(@() fclose(fid));
    for k=1:16
        if feof(fid), return; end
        fgetl(fid);
    end
    C = textscan(fid, '%s%s%s%s', 'Delimiter', ',', ...
        'ReturnOnError', false, 'MultipleDelimsAsOne', true, ...
        'Whitespace','', 'EndOfLine','\n');
    if isempty(C) || numel(C) < 4, return; end
    T = table(string(C{1}), string(C{2}), string(C{3}), string(C{4}));
end

function [rmse, mae, mape] = error_metrics(base_t, base_mov, ov_t, ov_v, rescale_flag, method)
    % Considera apenas sobreposição temporal
    in_rng = (ov_t >= min(base_t)) & (ov_t <= max(base_t));
    if nnz(in_rng) < 2
        rmse = NaN; mae = NaN; mape = NaN; return;
    end
    t_ref  = ov_t(in_rng);
    v_ov   = ov_v(in_rng);
    v_base = interp1(base_t, base_mov, t_ref, 'linear', 'extrap');

    if rescale_flag
        switch lower(method)
            case 'range'
                vmin = min(v_base,[],'omitnan'); vmax = max(v_base,[],'omitnan');
                if ~(isnan(vmin) || isnan(vmax)) && vmax > vmin
                    v_ov = rescale(v_ov, vmin, vmax);
                end
            case 'zscore'
                mu_b = mean(v_base,'omitnan'); sd_b = std(v_base,'omitnan');
                mu_o = mean(v_ov,'omitnan');   sd_o = std(v_ov,'omitnan');
                if sd_b > 0 && sd_o > 0
                    v_ov = ((v_ov - mu_o)/sd_o)*sd_b + mu_b;
                end
        end
    end

    e    = v_base - v_ov;
    rmse = sqrt(mean(e.^2,'omitnan'));
    mae  = mean(abs(e),'omitnan');
    epsm = 1e-6;
    mape = mean(abs(e) ./ max(abs(v_base), epsm), 'omitnan') * 100;
end

function [xmin, xmax, ymin, ymax] = xy_limits_union(t_base_mov, v_base_mov, t_ov, v_ov, v_ov_mov)
    xmin = min(t_base_mov); xmax = max(t_base_mov);
    ymin = min(v_base_mov,[],'omitnan'); ymax = max(v_base_mov,[],'omitnan');
    if ~isempty(t_ov)
        xmin = min([xmin; min(t_ov)]);
        xmax = max([xmax; max(t_ov)]);
    end
    if ~isempty(v_ov)
        ymin = min([ymin; min(v_ov,[],'omitnan')]);
        ymax = max([ymax; max(v_ov,[],'omitnan')]);
    end
    if ~isempty(v_ov_mov)
        ymin = min([ymin; min(v_ov_mov,[],'omitnan')]);
        ymax = max([ymax; max(v_ov_mov,[],'omitnan')]);
    end
    if ~isfinite(ymin) || ~isfinite(ymax) || ymin==ymax
        pad = 1;
        ymin = ymin - pad; ymax = ymax + pad;
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
