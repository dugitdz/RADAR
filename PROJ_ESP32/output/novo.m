clc; clear; close all;

radar_HR_path = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes.csv';
radar_RR_path = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\tes_rr.csv';

omni_HR_path = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\omni\Heart3.csv';
omni_RR_path = 'C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\omni\Breath3.csv';

[t1,RR_radar] = read_csv_radar(radar_RR_path);
[t2,HR_radar] = read_csv_radar(radar_HR_path);

[t3,RR_omni] = read_csv_omni(omni_RR_path);
[t4,HR_omni] = read_csv_omni(omni_HR_path);

t0_radar = min([t1(1),t2(1)]);
t0_omni = min([t3(1),t4(1)]);
 
mean_RR_radar = movmean(RR_radar,15);
mean_HR_radar = movmean(HR_radar,15);

mean_RR_omni = movmean(RR_omni,15);
mean_HR_omni = movmean(HR_omni,15);

figure;
plot((t1-t0_radar),mean_RR_radar,'red');
hold on;
plot(((double(t3-t0_omni))/1000),mean_RR_omni,'blue');
legend('Radar','Omni');

figure;
plot((t2-t0_radar),mean_HR_radar,'red');
hold on;
plot(((double(t4-t0_omni))/1000),mean_HR_omni,'blue');
legend('Radar','Omni');


%% ========================== CÁLCULO DO ERRO ==========================

% ---------- Preparar tempos ----------
t_RR_radar = (t1 - t0_radar);
t_RR_omni  = double(t3 - t0_omni) / 1000;

t_HR_radar = (t2 - t0_radar);
t_HR_omni  = double(t4 - t0_omni) / 1000;

% ---------- Criar grade comum ----------
dt_RR = min(median(diff(t_RR_radar)), median(diff(t_RR_omni)));
dt_HR = min(median(diff(t_HR_radar)), median(diff(t_HR_omni)));

tmin_RR = max([t_RR_radar(1), t_RR_omni(1)]);
tmax_RR = min([t_RR_radar(end), t_RR_omni(end)]);
t_grid_RR = (tmin_RR : dt_RR : tmax_RR).';

tmin_HR = max([t_HR_radar(1), t_HR_omni(1)]);
tmax_HR = min([t_HR_radar(end), t_HR_omni(end)]);
t_grid_HR = (tmin_HR : dt_HR : tmax_HR).';

% ---------- Interpolar ----------
RR_radar_grid = interp1(t_RR_radar, mean_RR_radar, t_grid_RR, 'linear');
RR_omni_grid  = interp1(t_RR_omni,  mean_RR_omni,  t_grid_RR, 'linear');

HR_radar_grid = interp1(t_HR_radar, mean_HR_radar, t_grid_HR, 'linear');
HR_omni_grid  = interp1(t_HR_omni,  mean_HR_omni,  t_grid_HR, 'linear');

% ---------- Erros ----------
err_RR = RR_radar_grid - RR_omni_grid;
err_HR = HR_radar_grid - HR_omni_grid;

MAE_RR  = mean(abs(err_RR));
RMSE_RR = sqrt(mean(err_RR.^2));

MAE_HR  = mean(abs(err_HR));
RMSE_HR = sqrt(mean(err_HR.^2));

% área entre curvas
area_RR = trapz(t_grid_RR, abs(err_RR));
area_HR = trapz(t_grid_HR, abs(err_HR));

% erro contínuo médio
mean_cont_RR = area_RR / (t_grid_RR(end) - t_grid_RR(1));
mean_cont_HR = area_HR / (t_grid_HR(end) - t_grid_HR(1));

% ---------- PRINT ----------
fprintf("\n========= ERROS RR =========\n");
fprintf("MAE: %.4f\n", MAE_RR);
fprintf("RMSE: %.4f\n", RMSE_RR);


fprintf("\n========= ERROS HR =========\n");
fprintf("MAE: %.4f\n", MAE_HR);
fprintf("RMSE: %.4f\n", RMSE_HR);

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
