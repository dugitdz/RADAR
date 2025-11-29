clc; clear; close all;

%% ========================== PATHS ==========================
radar_HR_path  = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes.csv';
radar_RR_path  = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes_rr.csv';
omni_HR_path   = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\omni\Heart4.csv';
omni_RR_path   = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\omni\Breath4.csv';
polar_HR_path  = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\POLARH2.txt';

%% ========================== LEITURA ==========================
% Polar (Referência)
[t_polar, HR_polar] = read_txt_polar(polar_HR_path);

% Radar Phase (Sinal a testar)
load('freq_results.mat'); % Carrega 't', 'mean_h_cont'
t_radar = t;
HR_radar = mean_h_cont;

% Remover duplicatas no tempo do Polar
[t_polar, ia] = unique(t_polar, 'stable');
HR_polar = HR_polar(ia);

%% ========================== ALINHAMENTO E INTERPOLAÇÃO (GRID COMUM) ==========================
% 1. Definir intervalo de tempo comum (Intersecção)
t_start = max(t_radar(1), t_polar(1));
t_end   = min(t_radar(end), t_polar(end));

% 2. Criar vetor de tempo regular (0.1s de passo) para cálculo de erro
dt_grid = 0.1; 
t_common = (t_start : dt_grid : t_end)';

% 3. Interpolar AMBOS os sinais para esse grid comum
HR_radar_grid = interp1(t_radar, HR_radar, t_common, 'linear');
HR_polar_grid = interp1(t_polar, HR_polar, t_common, 'linear');

%% ========================== FILTRAGEM DE INÍCIO (Opcional) ==========================
IGNORE_START_SEC = 5; % Ignora os primeiros 5 segundos para estabilização

mask_time = (t_common >= (t_start + IGNORE_START_SEC));
HR_ref  = HR_polar_grid(mask_time);
HR_test = HR_radar_grid(mask_time);
t_eval  = t_common(mask_time);

%% ========================== CÁLCULO DE ERRO E CORRELAÇÃO ==========================
error_vec = HR_test - HR_ref;

MAE  = mean(abs(error_vec), 'omitnan');
RMSE = sqrt(mean(error_vec.^2, 'omitnan'));
R    = corr(HR_test, HR_ref, 'rows', 'complete');

%% ========================== PLOT ÚNICO (TEMPORAL) ==========================
figure('Color','w','Position',[100 100 800 500]);

% Apenas o gráfico temporal
plot(t_polar, HR_polar, 'm.-', 'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(t_radar, HR_radar, 'g-', 'LineWidth', 1.5);

xlim([t_start t_end]);
title('Comparação Temporal: Radar Phase vs Polar');
xlabel('Tempo [s]'); ylabel('HR [bpm]');
legend('Polar (Referência)', 'Radar (Processado)');
grid on;

%% ========================== RESULTADOS ==========================

fprintf('MAE : %.4f bpm\n', MAE);
fprintf('RMSE: %.4f bpm\n', RMSE);
fprintf('Corr: %.4f\n', R);

%% ========================== FUNÇÕES AUXILIARES ==========================
function [t_sec, HR] = read_txt_polar(p)
    fid = fopen(p,'r');
    if fid == -1, error('Não foi possível abrir o arquivo Polar.'); end
    C = textscan(fid,'%s %f %*[^\n]','Delimiter',';','HeaderLines',1);
    fclose(fid);
    ts_str = C{1};
    HR = C{2};
    try
        t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
    catch
        t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    end
    t_sec = seconds(t_dt - t_dt(1));
end