phases_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv';

data   = readmatrix(phases_path,'CommentStyle', '#', 'Delimiter', ',');
tpuro  = data(:,1);
total  = data(:,2);
breath = data(:,3);
heart  = data(:,4);

N=128;
i=1;


t = (tpuro - tpuro(1)) / 1000;    % tempo em segundos
nt = length(t);

figure
plot(t,total,'black');
title('TOTAL');

figure
plot(t,breath,'b');
title('BREATH');

figure
plot(t,heart,'r');
title('HEART');
freq_parcial_h = zeros(size(t));
freq_parcial_br = zeros(size(t));
aux_h  = zeros(size(t));
aux_br = zeros(size(t));
while i*N<nt
com = 1+((i-2)*N/2);
if com<=0 com = 1; end
fim = i*N/2;
temp_novo = t(com:fim);
br_novo = breath(com:fim);
h_novo = heart(com:fim);

fs = 1/mean(diff(temp_novo));

amostra = (fim - com + 1)/2;

f = fft(h_novo);
fft_int_h = abs(f(1:amostra));
f = fft(br_novo);
fft_int_br = abs(f(1:amostra));

[~,maior] = max(fft_int_h);
freq_parcial_h(com:fim) = (maior-1)*fs/amostra;
[~,maior] = max(fft_int_br);
freq_parcial_br(com:fim) = (maior-1)*fs/amostra;

freq_h = (freq_parcial_h +  aux_h)/2;
freq_br = (freq_parcial_br + aux_br)/2;

aux_h = freq_parcial_h;
aux_br = freq_parcial_br;

i=i+1;
end




