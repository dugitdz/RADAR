clc; clear; close all;

%% ========================== PATHS ==========================
radar_HR_path  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes.csv';
radar_RR_path  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes_rr.csv';

omni_HR_path   = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\omni\Heart4.csv';
omni_RR_path   = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\omni\Breath4.csv';

polar_HR_path  = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH10.txt';

%% ========================== LEITURA DOS DADOS ==========================

% Radar
[t1,RR_radar] = read_csv_radar(radar_RR_path);
[t2,HR_radar] = read_csv_radar(radar_HR_path);

% Omni
[t3,RR_omni] = read_csv_omni(omni_RR_path);
[t4,HR_omni] = read_csv_omni(omni_HR_path);

% Polar (só HR)
[t5,HR_polar] = read_txt_polar(polar_HR_path);   % t5 em segundos desde o primeiro sample

%% ========================== ALINHO INICIAL ==========================

t0_radar = min([t1(1),t2(1)]);
t0_omni  = min([t3(1),t4(1)]);
% Polar já vem normalizado em segundos (começando em 0), não precisa de t0

%% ========================== SUAVIZAÇÃO ==========================

mean_RR_radar = movmean(RR_radar,15);
mean_HR_radar = movmean(HR_radar,15);

mean_RR_omni  = movmean(RR_omni,15);
mean_HR_omni  = movmean(HR_omni,15);

% Polar SEM movmean (como você pediu)
HR_polar_smooth = HR_polar;  % se quiser testar suavizado, é só trocar aqui

%% ========================== PLOTS ==========================

% RR: Radar x Omni
figure;
plot((t1-t0_radar), mean_RR_radar, 'red');
hold on;
plot((double(t3-t0_omni))/1000, mean_RR_omni, 'blue');
legend('Radar','Omni');
xlabel('Tempo [s]');
ylabel('RR [rpm]');
title('RR - Radar x Omni');

% HR: Radar x Omni x Polar
figure;
plot((t2-t0_radar), mean_HR_radar, 'red');
hold on;
plot((double(t4-t0_omni))/1000, mean_HR_omni, 'blue');
plot(t5, HR_polar_smooth, 'green');
legend('Radar','Omni','Polar');
xlabel('Tempo [s]');
ylabel('HR [bpm]');
title('HR - Radar x Omni x Polar');

%% ========================== CÁLCULO DO ERRO ==========================
% ---------- Preparar tempos ----------

t_RR_radar = (t1 - t0_radar);                    % s
t_RR_omni  = double(t3 - t0_omni) / 1000;        % s

t_HR_radar = (t2 - t0_radar);                    % s
t_HR_omni  = double(t4 - t0_omni) / 1000;        % s
t_HR_polar = t5;                                 % s (já normalizado na leitura)

%% ========================== ERROS RR (Radar x Omni) ==========================

dt_RR = min(median(diff(t_RR_radar)), median(diff(t_RR_omni)));

tmin_RR = max([t_RR_radar(1), t_RR_omni(1)]);
tmax_RR = min([t_RR_radar(end), t_RR_omni(end)]);
t_grid_RR = (tmin_RR : dt_RR : tmax_RR).';

RR_radar_grid = interp1(t_RR_radar, mean_RR_radar, t_grid_RR, 'linear');
RR_omni_grid  = interp1(t_RR_omni,  mean_RR_omni,  t_grid_RR, 'linear');

err_RR = RR_radar_grid - RR_omni_grid;

MAE_RR  = mean(abs(err_RR));
RMSE_RR = sqrt(mean(err_RR.^2));

area_RR     = trapz(t_grid_RR, abs(err_RR));
mean_cont_RR = area_RR / (t_grid_RR(end) - t_grid_RR(1));

%% ========================== ERROS HR (Radar x Omni) ==========================

dt_HR_RO = min(median(diff(t_HR_radar)), median(diff(t_HR_omni)));

tmin_HR_RO = max([t_HR_radar(1), t_HR_omni(1)]);
tmax_HR_RO = min([t_HR_radar(end), t_HR_omni(end)]);
t_grid_HR_RO = (tmin_HR_RO : dt_HR_RO : tmax_HR_RO).';

HR_radar_RO = interp1(t_HR_radar, mean_HR_radar, t_grid_HR_RO, 'linear');
HR_omni_RO  = interp1(t_HR_omni,  mean_HR_omni,  t_grid_HR_RO, 'linear');

err_HR_omni = HR_radar_RO - HR_omni_RO;

MAE_HR_omni  = mean(abs(err_HR_omni));
RMSE_HR_omni = sqrt(mean(err_HR_omni.^2));

area_HR_omni      = trapz(t_grid_HR_RO, abs(err_HR_omni));
mean_cont_HR_omni = area_HR_omni / (t_grid_HR_RO(end) - t_grid_HR_RO(1));

%% ========================== ERROS HR (Radar x Polar) ==========================

dt_HR_RP = min(median(diff(t_HR_radar)), median(diff(t_HR_polar)));

tmin_HR_RP = max([t_HR_radar(1), t_HR_polar(1)]);
tmax_HR_RP = min([t_HR_radar(end), t_HR_polar(end)]);
t_grid_HR_RP = (tmin_HR_RP : dt_HR_RP : tmax_HR_RP).';

HR_radar_RP = interp1(t_HR_radar, mean_HR_radar, t_grid_HR_RP, 'linear');
HR_polar_RP = interp1(t_HR_polar, HR_polar_smooth, t_grid_HR_RP, 'linear');

err_HR_polar = HR_radar_RP - HR_polar_RP;

MAE_HR_polar  = mean(abs(err_HR_polar));
RMSE_HR_polar = sqrt(mean(err_HR_polar.^2));

area_HR_polar      = trapz(t_grid_HR_RP, abs(err_HR_polar));
mean_cont_HR_polar = area_HR_polar / (t_grid_HR_RP(end) - t_grid_HR_RP(1));

%% ========================== PRINT ==========================

fprintf("\n========= ERROS RR (Radar x Omni) =========\n");
fprintf("MAE:  %.4f\n", MAE_RR);
fprintf("RMSE: %.4f\n", RMSE_RR);
fprintf("Erro contínuo médio (área/t): %.4f\n", mean_cont_RR);

fprintf("\n========= ERROS HR (Radar x Omni) =========\n");
fprintf("MAE:  %.4f\n", MAE_HR_omni);
fprintf("RMSE: %.4f\n", RMSE_HR_omni);
fprintf("Erro contínuo médio (área/t): %.4f\n", mean_cont_HR_omni);

fprintf("\n========= ERROS HR (Radar x Polar) =========\n");
fprintf("MAE:  %.4f\n", MAE_HR_polar);
fprintf("RMSE: %.4f\n", RMSE_HR_polar);
fprintf("Erro contínuo médio (área/t): %.4f\n", mean_cont_HR_polar);

%% ========================== FUNÇÕES ==========================

function [t, v] = read_csv_radar(p)
    data = readmatrix(p, 'CommentStyle','#', 'Delimiter', ',');
    t = data(:,1);
    v = data(:,2);
end

function [t, v] = read_csv_omni(p)
    data = readcell(p, 'CommentStyle','#', 'Delimiter', ',');
    v_raw = data(:,3);
    t_raw = data(:,6);

    if isnan(str2double(string(v_raw{1})))
        v_raw = v_raw(2:end);
        t_raw = t_raw(2:end);
    end

    v = cellfun(@(x) str2double(string(x)), v_raw);
    t = cellfun(@(x) str2double(erase(string(x),"'")), t_raw);
    t = int64(t);
end

function [t_sec, HR] = read_txt_polar(p)
    % Lê arquivo Polar .txt:
    % Phone timestamp;HR [bpm];HRV [ms];Breathing interval [rpm];
    fid = fopen(p,'r');
    if fid < 0
        error('Não foi possível abrir o arquivo Polar: %s', p);
    end

    % 1a coluna: string do timestamp
    % 2a coluna: HR (float/int)
    % resto da linha é ignorado
    C = textscan(fid, '%s %f %*[^\n]', 'Delimiter',';', 'HeaderLines',1);
    fclose(fid);

    ts_str = C{1};
    HR     = C{2};

    % Converte para datetime (formato 2025-11-17T17:26:16.058)
    t_dt = datetime(ts_str, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');

    % Converte para segundos desde o primeiro sample
    t_sec = seconds(t_dt - t_dt(1));
end
