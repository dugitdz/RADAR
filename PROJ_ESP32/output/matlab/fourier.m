phases_path = 'C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv';

data   = readmatrix(phases_path,'CommentStyle', '#', 'Delimiter', ',');
tpuro  = data(:,1);
total  = data(:,2);
breath = data(:,3);
heart  = data(:,4);

N = 128;
i = 1;

t  = (tpuro - tpuro(1)) / 1000;
nt = length(t);

figure
plot(t, total, 'black');
title('TOTAL');

figure
plot(t, breath, 'b');
title('BREATH');

figure
plot(t, heart, 'r');
title('HEART');

freq_parcial_h  = zeros(size(t));
freq_parcial_br = zeros(size(t));
aux_h  = zeros(size(t));
aux_br = zeros(size(t));
freq_h  = zeros(size(t));
freq_br = zeros(size(t));

while (i*N/2) <= nt

    com = 1 + ((i-2)*N/2);
    if com <= 0
        com = 1;
    end
    fim = i*N/2;

    if fim > nt
        fim = nt;
    end

    temp_novo = t(com:fim);
    br_novo   = breath(com:fim);
    h_novo    = heart(com:fim);

    if numel(temp_novo) < 2
        break;
    end

    dt_local = diff(temp_novo);
    dt_medio = mean(dt_local);
    dt_max   = max(dt_local);

    if dt_medio > 0.1 || dt_max > 0.5
        freq_parcial_h(com:fim)  = aux_h(com:fim);
        freq_parcial_br(com:fim) = aux_br(com:fim);
    else
        fs = 1/dt_medio;

        L = fim - com + 1;
        amostra = floor(L/2);

        f = fft(h_novo);
        fft_int_h = abs(f(1:amostra));
        if amostra >= 1
            fft_int_h(1) = 0;
        end
        [~, maior] = max(fft_int_h);
        freq_parcial_h(com:fim) = (maior-1)*fs / L;

        f = fft(br_novo);
        fft_int_br = abs(f(1:amostra));
        if amostra >= 1
            fft_int_br(1) = 0;
        end
        [~, maior] = max(fft_int_br);
        freq_parcial_br(com:fim) = (maior-1)*fs / L;
    end

    if i == 1
        freq_h  = freq_parcial_h;
        freq_br = freq_parcial_br;
    else
        freq_h  = freq_parcial_h;
        freq_br = freq_parcial_br;

        idx_h = (aux_h ~= 0) & (freq_parcial_h ~= 0);
        freq_h(idx_h) = (freq_parcial_h(idx_h) + aux_h(idx_h)) / 2;

        idx_br = (aux_br ~= 0) & (freq_parcial_br ~= 0);
        freq_br(idx_br) = (freq_parcial_br(idx_br) + aux_br(idx_br)) / 2;
    end 

    aux_h  = freq_h;
    aux_br = freq_br;

    i = i + 1;
end

idx = find(freq_h ~= 0);
if ~isempty(idx)
    first = idx(1);
    last  = idx(end);
    if first > 1
        freq_h(1:first-1) = freq_h(first);
    end
    if last < nt
        freq_h(last+1:end) = freq_h(last);
    end
    i = first + 1;
    while i < last
        if freq_h(i) == 0 && freq_h(i-1) ~= 0
            j = i + 1;
            while j <= last && freq_h(j) == 0
                j = j + 1;
            end
            if j <= last
                for k = i:j-1
                    freq_h(k) = freq_h(i-1) + (freq_h(j) - freq_h(i-1)) * (k - (i-1)) / (j - (i-1));
                end
                i = j;
                continue;
            else
                break;
            end
        end
        i = i + 1;
    end
end

idx = find(freq_br ~= 0);
if ~isempty(idx)
    first = idx(1);
    last  = idx(end);
    if first > 1
        freq_br(1:first-1) = freq_br(first);
    end
    if last < nt
        freq_br(last+1:end) = freq_br(last);
    end
    i = first + 1;
    while i < last
        if freq_br(i) == 0 && freq_br(i-1) ~= 0
            j = i + 1;
            while j <= last && freq_br(j) == 0
                j = j + 1;
            end
            if j <= last
                for k = i:j-1
                    freq_br(k) = freq_br(i-1) + (freq_br(j) - freq_br(i-1)) * (k - (i-1)) / (j - (i-1));
                end
                i = j;
                continue;
            else
                break;
            end
        end
        i = i + 1;
    end
end

freq_br = 60*freq_br;
freq_h  = 60*freq_h;

mean_h  = movmean(freq_h,15);
mean_br = movmean(freq_br,15);

figure
plot(t,freq_br,'r');
hold on;
plot(t,mean_br,'b');

figure
plot(t,freq_h,'r');
hold on;
plot(t,mean_h,'b');
save('freq_results.mat','freq_h','freq_br','mean_h','mean_br','t');
