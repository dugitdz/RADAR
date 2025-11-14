clc; clear; close all;

radar_HR_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes.csv';
radar_RR_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes_rr.csv';

omni_HR_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\omni\Heart1.csv';
omni_RR_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\omni\Breath1.csv';

[t1,RR_radar] = read_csv_radar(radar_RR_path);
[t2,HR_radar] = read_csv_radar(radar_HR_path);

[t3,RR_omni] = read_csv_omni(omni_RR_path);
[t4,HR_omni] = read_csv_omni(omni_HR_path);

t0_radar = min([t1(1),t2(1)]);
t0_omni = min([t3(1),t4(1)]);
 
plot((t1-t0_radar),RR_radar,'red');
hold on;
plot(((double(t3-t0_omni))/1000),RR_omni);



function [t, v] = read_csv_radar(p)
    if exist(p,'file') ~= 2
        error('Arquivo base não encontrado: %s', p);
    end
    data = readmatrix(p, 'CommentStyle','#', 'Delimiter', ',');
    if isempty(data) || size(data,2) < 2
        error('Base inválida (precisa "dt_seconds,value"): %s', p);
    end
    t = data(:,1);
    v = data(:,2);
 end

 function [t, v] = read_csv_omni(p)
    if exist(p,'file') ~= 2
        error('Arquivo base não encontrado: %s', p);
    end
    data = readcell(p, 'CommentStyle','#', 'Delimiter', ',');
   
    v = data(:,3);
    t_sujo = data(:,6);
    t = int64(cellfun(@(x) str2double(erase(x,"'")), t_sujo));
 end
