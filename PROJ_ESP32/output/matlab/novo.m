clc; clear; close all;

%% ===================== PARÂMETROS GERAIS ===================== %%

% banda fisiológica do coração (para FFT / detecção de pico nas janelas)
f_low  = 0.8;   % Hz
f_high = 3.0;   % Hz

% X fixo -> h = X / frex
X = 0.6;

% caminhos base
base_path = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\';

% ==================== 3 DUPLAS (PHASE x POLAR) ==================== %
pairs = {
    struct('phase',  'phase.csv',      'polar', 'POLARH10.txt' );
    struct('phase',  'phases.csv',     'polar', 'POLARH2.txt'  );
    struct('phase',  'phases_raw.csv', 'polar', 'POLARH3.txt'  );
};

nPairs = numel(pairs);

% ==================== GRID DE PARÂMETROS ==================== %
N_list         = 8:8:256;        % tamanho da janela FFT
overlap_list   = 0.5:0.05:0.95;   % sobreposição
win_raw_list   = 32:32:128;        % janelas de movmean (sinal bruto)
win_final_list = 32:32:128;        % janelas de movmean (sinal final)

% normalização do RMSE pra ficar na mesma ordem de (1 - corr)
rmse_scale = 10;   % bpm

best_score_global = Inf;
best_cfg = struct('N',[],'overlap',[],'win_raw',[],'win_final',[], ...
                  'rmse_pairs',[],'corr_pairs',[]);

%% ==================== LOOP DO GRID ==================== %%

for iN = 1:length(N_list)
    N = N_list(iN);

    for io = 1:length(overlap_list)
        overlap = overlap_list(io);

        for iwr = 1:length(win_raw_list)
            win_smooth_raw = win_raw_list(iwr);

            for iwf = 1:length(win_final_list)
                win_smooth_final = win_final_list(iwf);

                rmse_list  = nan(1,nPairs);
                corr_list  = nan(1,nPairs);
                score_list = nan(1,nPairs);

                fprintf('---- Testando N=%3d | over=%.2f | win_raw=%3d | win_final=%3d ----\n', ...
                        N, overlap, win_smooth_raw, win_smooth_final);

                for p = 1:nPairs
                    phase_path = fullfile(base_path, pairs{p}.phase);
                    polar_path = fullfile(base_path, pairs{p}.polar);

                    [rmse_p, corr_p] = compute_error_for_pair( ...
                        phase_path, polar_path, X, ...
                        N, overlap, win_smooth_raw, win_smooth_final, ...
                        f_low, f_high);

                    rmse_list(p) = rmse_p;
                    corr_list(p) = corr_p;

                    fprintf('   Par %d: RMSE = %.4f | Corr = %.4f\n', ...
                            p, rmse_p, corr_p);

                    if isnan(rmse_p) || isnan(corr_p)
                        score_list(p) = NaN;
                    else
                        rmse_norm = rmse_p / rmse_scale;
                        % SCORE balanceado: RMSE_norm + (1 - corr)
                        score_list(p) = rmse_norm + (1 - corr_p);
                    end
                end

                % ignora configs que falharam em algum par
                if any(isnan(score_list))
                    continue;
                end

                score_global = mean(score_list);

                if score_global < best_score_global
                    best_score_global = score_global;
                    best_cfg.N          = N;
                    best_cfg.overlap    = overlap;
                    best_cfg.win_raw    = win_smooth_raw;
                    best_cfg.win_final  = win_smooth_final;
                    best_cfg.rmse_pairs = rmse_list;
                    best_cfg.corr_pairs = corr_list;

                    fprintf('\n>>> Novo MELHOR!  SCORE = %.4f\n', score_global);
                    fprintf('    N=%3d | overlap=%.2f | win_raw=%3d | win_final=%3d\n\n', ...
                            N, overlap, win_smooth_raw, win_smooth_final);
                end

            end
        end
    end
end

%% ==================== RESULTADO FINAL ==================== %%

fprintf('\n================= RESULTADO FINAL DO GRID =================\n');
fprintf('Melhor configuração:\n');
fprintf('  N         = %d\n',   best_cfg.N);
fprintf('  overlap   = %.2f\n', best_cfg.overlap);
fprintf('  win_raw   = %d (movmean sinal bruto)\n',  best_cfg.win_raw);
fprintf('  win_final = %d (movmean sinal final)\n', best_cfg.win_final);
fprintf('  SCORE mín = %.4f  (média de [RMSE/%.1f + (1-corr)] nos 3 pares)\n\n', ...
        best_score_global, rmse_scale);

for p = 1:nPairs
    fprintf('Par %d (%s x %s): RMSE = %.4f | Corr = %.4f\n', ...
        p, pairs{p}.phase, pairs{p}.polar, ...
        best_cfg.rmse_pairs(p), best_cfg.corr_pairs(p));
end


%% ==================== FUNÇÃO: ERRO PARA 1 PAR (PHASE x POLAR) ==================== %%
function [RMSE_HR_PP, Corr_HR_PP] = compute_error_for_pair( ...
    phase_path, polar_path, X, ...
    N, overlap, win_smooth_raw, win_smooth_final, ...
    f_low, f_high)

    % -------- LEITURA DO CSV DE PHASE --------
    data  = readmatrix(phase_path);
    tpuro = data(:,1);
    heart = data(:,4);

    % tempo em segundos
    t = (tpuro - tpuro(1)) / 1000;
    dt = mean(diff(t));
    srate = 1/dt;
    nt = numel(t);

    % -------- DETECTAR QUEBRAS dt > 0.5 s E SEGMENTOS --------
    dt_all   = diff(t);
    break_ix = find(dt_all > 0.5);   % índices onde há "buraco" no tempo

    if isempty(break_ix)
        seg_starts = 1;
        seg_ends   = nt;
    else
        seg_starts = [1; break_ix + 1];
        seg_ends   = [break_ix; nt];
    end
    nSeg = numel(seg_starts);

    % -------- AUTO FREX PELO ESPECTRO DO HEART --------
    Nfft  = 2*nt - 1;
    f_fft = linspace(0, srate, Nfft);
    H     = fft(heart, Nfft);
    magH  = abs(H);

    [~, imax] = max(magH);
    frex = f_fft(imax);

    if frex <= 0
        RMSE_HR_PP = NaN;
        Corr_HR_PP = NaN;
        return;
    end

    % -------- DEFINIR h = X/frex --------
    h = X / frex;         % em segundos
    h = max(h, 0.1);      % limites de segurança
    h = min(h, 5.0);

    % -------- WAVELET & CONVOLUÇÃO --------
    idx   = (1:nt) - ceil(nt/2);
    t_aux = idx(:) * dt;

    gaus = exp((-4*log(2)*t_aux.^2)/h^2) .* exp(1i*2*pi*frex*t_aux);

    nGaus    = length(gaus);
    nHeart   = length(heart);
    nConv    = nGaus + nHeart - 1;

    heartX = fft(heart, nConv);
    gausX  = fft(gaus,  nConv);
    gausX  = gausX ./ max(gausX);

    convX     = heartX .* gausX;
    convs_all = ifft(convX);

    % recorte "same": mesmo tamanho e alinhado ao heart
    startIdx = floor(nGaus/2) + 1;
    endIdx   = startIdx + nt - 1;
    convs    = convs_all(startIdx:endIdx);

    % -------- FFT POR JANELAS EM "convs" --------
    hop     = max(1, round(N*(1-overlap)));
    halfN   = N/2;
    rel_idx = -halfN : (halfN-1);

    freq_weighted_sum = zeros(nt,1);
    weight_sum        = zeros(nt,1);

    for s = 1:nSeg
        seg_start = seg_starts(s);
        seg_end   = seg_ends(s);

        center_list = seg_start : hop : seg_end;

        for center_theoretical = center_list
            idx_theo = center_theoretical + rel_idx;
            idx_win  = idx_theo(idx_theo >= seg_start & idx_theo <= seg_end);

            if numel(idx_win) < 4
                continue;
            end

            t_seg = t(idx_win);
            x_seg = real(convs(idx_win));

            dt_seg = diff(t_seg);
            if any(dt_seg <= 0) || any(isnan(dt_seg))
                continue;
            end

            fs_seg = 1/mean(dt_seg);

            x_seg = x_seg - mean(x_seg);

            Xw        = fft(x_seg);
            nFFT_win  = numel(x_seg);
            halfFFT   = floor(nFFT_win/2);
            f_ax      = (0:halfFFT-1) * fs_seg / nFFT_win;
            magX_win  = abs(Xw(1:halfFFT));

            mask = (f_ax >= f_low) & (f_ax <= f_high);
            if ~any(mask)
                continue;
            end

            [~, imax_win] = max(magX_win(mask));
            f_band = f_ax(mask);
            f_est  = f_band(imax_win);
            bpm    = f_est * 60;

            % acumula contribuição ponderada desta janela
            for ii = idx_win
                dist = abs(ii - center_theoretical);
                w = (halfN + 1) - dist;
                if w < 0, w = 0; end
                freq_weighted_sum(ii) = freq_weighted_sum(ii) + bpm * w;
                weight_sum(ii)        = weight_sum(ii)        + w;
            end
        end
    end

    % -------- HR POR AMOSTRA (bruto) --------
    freq_h_raw = nan(nt,1);
    valid = weight_sum > 0;
    freq_h_raw(valid) = freq_weighted_sum(valid) ./ weight_sum(valid);

    % -------- SUAVIZAÇÃO DENTRO DE CADA SEGMENTO (RAW) --------
    freq_h_smooth = nan(nt,1);
    for s = 1:nSeg
        idx_seg = seg_starts(s):seg_ends(s);
        seg_raw = freq_h_raw(idx_seg);
        freq_h_smooth(idx_seg) = movmean(seg_raw, win_smooth_raw, 'omitnan');
    end

    % -------- INTERPOLAÇÃO PARA FECHAR BURACOS --------
    idx_good = ~isnan(freq_h_smooth);
    if any(idx_good) && any(~idx_good)
        freq_h_smooth(~idx_good) = interp1(t(idx_good), freq_h_smooth(idx_good), ...
                                           t(~idx_good), 'linear', 'extrap');
    end

    % -------- SEGUNDA SUAVIZAÇÃO (FINAL, POR SEGMENTO) --------
    freq_h_final = nan(nt,1);
    for s = 1:nSeg
        idx_seg     = seg_starts(s):seg_ends(s);
        seg_smooth  = freq_h_smooth(idx_seg);
        freq_h_final(idx_seg) = movmean(seg_smooth, win_smooth_final, 'omitnan');
    end

    t_phase        = t;
    HR_radar_phase = freq_h_final;

    % -------- LER POLAR --------
    [t_polar, HR_polar] = read_txt_polar(polar_path);
    [t_polar, ia] = unique(t_polar,'stable');
    HR_polar = HR_polar(ia);

    % -------- ALINHAR E COMPARAR --------
    dt_HR_PP   = min(median(diff(t_phase)), median(diff(t_polar)));
    tmin_HR_PP = max([t_phase(1), t_polar(1)]);
    tmax_HR_PP = min([t_phase(end), t_polar(end)]);

    if tmax_HR_PP <= tmin_HR_PP
        RMSE_HR_PP = NaN;
        Corr_HR_PP = NaN;
        return;
    end

    t_grid_HR_PP = (tmin_HR_PP : dt_HR_PP : tmax_HR_PP).';

    HR_phase_PP = interp1(t_phase, HR_radar_phase, t_grid_HR_PP,'linear');
    HR_polar_PP = interp1(t_polar, HR_polar,      t_grid_HR_PP,'linear');

    % descartar início (lock-in), se tiver muitos pontos
    if numel(HR_phase_PP) > 250
        HR_phase_PP = HR_phase_PP(251:end);
        HR_polar_PP = HR_polar_PP(251:end);
    end

    err = HR_phase_PP - HR_polar_PP;

    RMSE_HR_PP = sqrt(mean(err.^2));
    Corr_HR_PP = corr(HR_phase_PP, HR_polar_PP, 'rows','complete');
end

%% ==================== FUNÇÃO: LEITURA POLAR ==================== %%
function [t_sec, HR] = read_txt_polar(p)
    fid = fopen(p,'r');
    C = textscan(fid,'%s %f %*[^\n]','Delimiter',';','HeaderLines',1);
    fclose(fid);

    ts_str = C{1};
    HR     = C{2};

    t_dt  = datetime(ts_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
    t_sec = seconds(t_dt - t_dt(1));
end
