clc; clear; close all;

%% ========================== PATHS ==========================
radar_HR_path  = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes.csv';
radar_RR_path  = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes_rr.csv';

omni_HR_path   = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\omni\Heart4.csv';
omni_RR_path   = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\omni\Breath4.csv';

polar_HR_path  = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\POLARH10.txt';

%% ========================== LEITURA DOS DADOS ==========================

% Radar (pronto)
[t1,RR_radar] = read_csv_radar(radar_RR_path);
[t2,HR_radar] = read_csv_radar(radar_HR_path);

% Omni
[t3,RR_omni] = read_csv_omni(omni_RR_path);
[t4,HR_omni] = read_csv_omni(omni_HR_path);

% Polar (só HR)
[t5,HR_polar] = read_txt_polar(polar_HR_path);

% Radar_phase vindo do .mat
load('freq_results.mat');     % carrega t, mean_h, mean_br
t_phase             = t;       % tempo em segundos
RR_radar_phase      = mean_br; % já suavizado no outro código
HR_radar_phase      = mean_h;

%% ========================== SUAVIZAÇÃO ==========================

mean_RR_radar = movmean(RR_radar,15);
mean_HR_radar = movmean(HR_radar,15);

mean_RR_omni  = movmean(RR_omni,15);
mean_HR_omni  = movmean(HR_omni,15);

% Polar SEM movmean — já definido
% Radar_phase SEM movmean — já definido

% ======== REMOVER TIMESTAMPS DUPLICADOS DO POLAR ==========
t_HR_polar = t5;
[t_HR_polar, ia] = unique(t_HR_polar,'stable');
HR_polar = HR_polar(ia);

%% ========================== ALINHO INICIAL ==========================

t0_radar = min([t1(1),t2(1)]);
t0_omni  = min([t3(1),t4(1)]);

%% ========================== PLOTS ==========================

% RR: Radar x Omni x Radar_phase
figure;
plot((t1 - t0_radar), mean_RR_radar, 'red');
hold on;
plot((double(t3 - t0_omni)/1000), mean_RR_omni, 'blue');
plot(t_phase, RR_radar_phase, 'green');
legend('Radar','Omni','Radar phase');
xlabel('Tempo [s]');
ylabel('RR [rpm]');
title('RR - Radar x Omni x Radar phase');

% HR: Radar x Omni x Polar x Radar_phase
figure;
plot((t2 - t0_radar), mean_HR_radar, 'red');
hold on;
plot((double(t4 - t0_omni)/1000), mean_HR_omni, 'blue');
plot(t_HR_polar, HR_polar, 'magenta');
plot(t_phase, HR_radar_phase, 'green');
legend('Radar','Omni','Polar','Radar phase');
xlabel('Tempo [s]');
ylabel('HR [bpm]');
title('HR - Radar x Omni x Polar x Radar phase');

%% ========================== CÁLCULO DO ERRO ==========================

% --- tempos RR ---
t_RR_radar       = (t1 - t0_radar);
t_RR_omni        = double(t3 - t0_omni)/1000;
t_RR_radar_phase = t_phase;

% --- tempos HR ---
t_HR_radar       = (t2 - t0_radar);
t_HR_omni        = double(t4 - t0_omni)/1000;
t_HR_polar       = t_HR_polar;
t_HR_radar_phase = t_phase;

%% ========================== RR: Radar x Omni ==========================

dt_RR = min(median(diff(t_RR_radar)), median(diff(t_RR_omni)));
tmin_RR = max([t_RR_radar(1), t_RR_omni(1)]);
tmax_RR = min([t_RR_radar(end), t_RR_omni(end)]);
t_grid_RR = (tmin_RR : dt_RR : tmax_RR).';

RR_radar_g = interp1(t_RR_radar, mean_RR_radar, t_grid_RR,'linear');
RR_omni_g  = interp1(t_RR_omni,  mean_RR_omni,  t_grid_RR,'linear');

err_RR = RR_radar_g - RR_omni_g;

MAE_RR  = mean(abs(err_RR));
RMSE_RR = sqrt(mean(err_RR.^2));

%% ========================== RR: Radar x Radar_phase ==========================

dt_RR_RP = min(median(diff(t_RR_radar)), median(diff(t_RR_radar_phase)));

tmin_RR_RP = max([t_RR_radar(1), t_RR_radar_phase(1)]);
tmax_RR_RP = min([t_RR_radar(end), t_RR_radar_phase(end)]);
t_grid_RR_RP = (tmin_RR_RP : dt_RR_RP : tmax_RR_RP).';

RR_radar_g_RP       = interp1(t_RR_radar, mean_RR_radar, t_grid_RR_RP,'linear');
RR_radar_phase_g_RP = interp1(t_RR_radar_phase, RR_radar_phase, t_grid_RR_RP,'linear');

err_RR_RP  = RR_radar_g_RP - RR_radar_phase_g_RP;
MAE_RR_RP  = mean(abs(err_RR_RP));
RMSE_RR_RP = sqrt(mean(err_RR_RP.^2));

%% ========================== RR: Radar_phase x Omni ==========================

dt_RR_PO = min(median(diff(t_RR_radar_phase)), median(diff(t_RR_omni)));

tmin_RR_PO = max([t_RR_radar_phase(1), t_RR_omni(1)]);
tmax_RR_PO = min([t_RR_radar_phase(end), t_RR_omni(end)]);
t_grid_RR_PO = (tmin_RR_PO : dt_RR_PO : tmax_RR_PO).';

RR_phase_PO = interp1(t_RR_radar_phase, RR_radar_phase, t_grid_RR_PO,'linear');
RR_omni_PO  = interp1(t_RR_omni,        mean_RR_omni,  t_grid_RR_PO,'linear');

err_RR_PO  = RR_phase_PO - RR_omni_PO;
MAE_RR_PO  = mean(abs(err_RR_PO));
RMSE_RR_PO = sqrt(mean(err_RR_PO.^2));

%% ========================== HR: Radar x Omni ==========================

dt_HR_RO = min(median(diff(t_HR_radar)), median(diff(t_HR_omni)));
tmin_HR_RO = max([t_HR_radar(1), t_HR_omni(1)]);
tmax_HR_RO = min([t_HR_radar(end), t_HR_omni(end)]);
t_grid_HR_RO = (tmin_HR_RO : dt_HR_RO : tmax_HR_RO).';

HR_radar_g = interp1(t_HR_radar, mean_HR_radar, t_grid_HR_RO,'linear');
HR_omni_g  = interp1(t_HR_omni,  mean_HR_omni,  t_grid_HR_RO,'linear');

err_HR_omni = HR_radar_g - HR_omni_g;

MAE_HR_omni  = mean(abs(err_HR_omni));
RMSE_HR_omni = sqrt(mean(err_HR_omni.^2));

%% ========================== HR: Radar x Polar ==========================

dt_HR_RP = min(median(diff(t_HR_radar)), median(diff(t_HR_polar)));
tmin_HR_RP = max([t_HR_radar(1), t_HR_polar(1)]);
tmax_HR_RP = min([t_HR_radar(end), t_HR_polar(end)]);
t_grid_HR_RP = (tmin_HR_RP : dt_HR_RP : tmax_HR_RP).';

HR_radar_RP = interp1(t_HR_radar, mean_HR_radar, t_grid_HR_RP,'linear');
HR_polar_RP = interp1(t_HR_polar, HR_polar,     t_grid_HR_RP,'linear');

err_HR_polar = HR_radar_RP - HR_polar_RP;

MAE_HR_polar  = mean(abs(err_HR_polar));
RMSE_HR_polar = sqrt(mean(err_HR_polar.^2));

%% ========================== HR: Radar x Radar_phase ==========================

dt_HR_RPH = min(median(diff(t_HR_radar)), median(diff(t_HR_radar_phase)));
tmin_HR_RPH = max([t_HR_radar(1), t_HR_radar_phase(1)]);
tmax_HR_RPH = min([t_HR_radar(end), t_HR_radar_phase(end)]);
t_grid_HR_RPH = (tmin_HR_RPH : dt_HR_RPH : tmax_HR_RPH).';

HR_radar_g2     = interp1(t_HR_radar,       mean_HR_radar,       t_grid_HR_RPH,'linear');
HR_phase_g2     = interp1(t_HR_radar_phase, HR_radar_phase,      t_grid_HR_RPH,'linear');

err_HR_phase  = HR_radar_g2 - HR_phase_g2;

MAE_HR_phase  = mean(abs(err_HR_phase));
RMSE_HR_phase = sqrt(mean(err_HR_phase.^2));

%% ========================== HR: Radar_phase x Omni ==========================

dt_HR_PO = min(median(diff(t_HR_radar_phase)), median(diff(t_HR_omni)));

tmin_HR_PO = max([t_HR_radar_phase(1), t_HR_omni(1)]);
tmax_HR_PO = min([t_HR_radar_phase(end), t_HR_omni(end)]);
t_grid_HR_PO = (tmin_HR_PO : dt_HR_PO : tmax_HR_PO).';

HR_phase_PO = interp1(t_HR_radar_phase, HR_radar_phase, t_grid_HR_PO,'linear');
HR_omni_PO  = interp1(t_HR_omni,        mean_HR_omni,   t_grid_HR_PO,'linear');

err_HR_PO = HR_phase_PO - HR_omni_PO;

MAE_HR_PO  = mean(abs(err_HR_PO));
RMSE_HR_PO = sqrt(mean(err_HR_PO.^2));

%% ========================== HR: Radar_phase x Polar ==========================

dt_HR_PP = min(median(diff(t_HR_radar_phase)), median(diff(t_HR_polar)));

tmin_HR_PP = max([t_HR_radar_phase(1), t_HR_polar(1)]);
tmax_HR_PP = min([t_HR_radar_phase(end), t_HR_polar(end)]);
t_grid_HR_PP = (tmin_HR_PP : dt_HR_PP : tmax_HR_PP).';

HR_phase_PP = interp1(t_HR_radar_phase, HR_radar_phase, t_grid_HR_PP,'linear');
HR_polar_PP = interp1(t_HR_polar,       HR_polar,        t_grid_HR_PP,'linear');

err_HR_PP = HR_phase_PP - HR_polar_PP;

MAE_HR_PP  = mean(abs(err_HR_PP));
RMSE_HR_PP = sqrt(mean(err_HR_PP.^2));

%% ========================== PRINT ==========================

fprintf("\n========= ERROS RR (Radar x Omni) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_RR, RMSE_RR);

fprintf("\n========= ERROS RR (Radar x Radar phase) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_RR_RP, RMSE_RR_RP);

fprintf("\n========= ERROS RR (Radar phase x Omni) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_RR_PO, RMSE_RR_PO);

fprintf("\n========= ERROS HR (Radar x Omni) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_HR_omni, RMSE_HR_omni);

fprintf("\n========= ERROS HR (Radar x Polar) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_HR_polar, RMSE_HR_polar);

fprintf("\n========= ERROS HR (Radar x Radar phase) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_HR_phase, RMSE_HR_phase);

fprintf("\n========= ERROS HR (Radar_phase x Omni) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_HR_PO, RMSE_HR_PO);

fprintf("\n========= ERROS HR (Radar_phase x Polar) =========\n");
fprintf("MAE: %.4f\nRMSE: %.4f\n", MAE_HR_PP, RMSE_HR_PP);

%% ========================== FUNÇÕES ==========================

function [t, v] = read_csv_radar(p)
    data = readmatrix(p,'CommentStyle','#','Delimiter',',');
    t = data(:,1); v = data(:,2);
end

function [t, v] = read_csv_omni(p)
    data = readcell(p,'CommentStyle','#','Delimiter',',');
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
    fid = fopen(p,'r');
    C = textscan(fid,'%s %f %*[^\n]','Delimiter',';','HeaderLines',1);
    fclose(fid);

    ts_str = C{1};
    HR = C{2};

    t_dt = datetime(ts_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');
    t_sec = seconds(t_dt - t_dt(1));
end
