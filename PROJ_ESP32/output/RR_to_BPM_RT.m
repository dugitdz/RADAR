clearvars; close all; clc;

% HEALTH RADAR v8
% Config fixa: Fs=50/3, DWT db5 L3 D3, winN=32, hop=16, mm=10s
% Confianca mantida (pico/soma, pico/mediana, pico/2o_pico). Sem jump penalty.
% NFFT=512 para resolucao espectral com janela curta.

BASE = "C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\";
DATASETS = [
    struct("name","H2",    "radar", BASE+"radar7.csv",  "polar", BASE+"polar7.txt")
    struct("name","H3",    "radar", BASE+"radar8.csv",  "polar", BASE+"polar8.txt")
    struct("name","H10",   "radar", BASE+"radar9.csv",  "polar", BASE+"polar9.txt")
    struct("name","TESTE", "radar", BASE+"radar10.csv", "polar", BASE+"polar10.txt")
];

COL_T_MS=1; COL_PHASE=2;
FS_NATIVE=50.0;
FS_TARGET=FS_NATIVE/3;   % 16.6667 Hz
DWT_WAVE='db5';
DWT_LEVEL=3;
KEEPD=3;                  % D3: 1.04-2.08 Hz = 62.5-125 bpm
winN=32;                  % 1.92 s
hopN=16;                  % 0.96 s
NFFT_BIG=512;             % bin = 0.033 Hz ~ 2 bpm
FMIN_HZ=0.9;
FMAX_HZ=2.0;
CONF_THR=0.105;          % threshold de confianca composta
FINAL_MOVMEAN_SEC=9;
gap_thr=0.5;
WRAP_ON=0;
MAX_LAG_SEC=10.0;
FILL_GAPS_FOR_PLOT=1; MAX_GAP_SEC_PLOT=4.0;

% Info
fprintf('\n======================================================================\n');
fprintf('  HEALTH RADAR v8 — Config fixa (grid result)\n');
fprintf('  Fs=%.4f Hz (50/3) | DWT %s L%d | D%d: %.4f-%.4f Hz (%.0f-%.0f bpm)\n',...
    FS_TARGET,DWT_WAVE,DWT_LEVEL,KEEPD,...
    FS_TARGET/2^(KEEPD+1),FS_TARGET/2^KEEPD,...
    FS_TARGET/2^(KEEPD+1)*60,FS_TARGET/2^KEEPD*60);
fprintf('  winN=%d (%.2fs) hop=%d (%.2fs) NFFT=%d (bin=%.4f Hz = %.1f bpm)\n',...
    winN,winN/FS_TARGET,hopN,hopN/FS_TARGET,NFFT_BIG,...
    FS_TARGET/NFFT_BIG,FS_TARGET/NFFT_BIG*60);
fprintf('  HR search: [%.2f-%.2f] Hz = [%d-%d] bpm | movmean=%ds\n',...
    FMIN_HZ,FMAX_HZ,round(FMIN_HZ*60),round(FMAX_HZ*60),FINAL_MOVMEAN_SEC);
fprintf('  Confianca: pico/soma + pico/mediana + pico/2o_pico | thr=%.2f\n',CONF_THR);
fprintf('======================================================================\n');

% Figura
fig=figure('Color','w','Position',[50 50 1100 850]);
tlay=tiledlayout(numel(DATASETS),1,'TileSpacing','tight','Padding','compact');
ax=gobjects(numel(DATASETS),1);
hPL=[];hNCL=[];

for d=1:numel(DATASETS)
    ds=DATASETS(d);
    rp=resolve_any_path(ds.radar);rl=get_file_stem(rp);
    pp=resolve_any_path(ds.polar);
    [tP,hP]=read_txt_polar_flex(pp);
    if isempty(tP),fprintf('\n%s | Polar vazio.\n',rl);continue;end
    nP=numel(tP);
    A_raw=[];
    try A_raw=readmatrix(rp);t_dur=(A_raw(end,COL_T_MS)-A_raw(1,COL_T_MS))/1000;
    catch,t_dur=NaN;end
    fprintf('\n-> [%s] Radar=%.1fs Polar=%.1fs (N=%d)\n',ds.name,t_dur,tP(end),nP);

    % == NON-CAUSAL ==
    [hNC,okNC]=algo_noncausal(rp,tP,COL_T_MS,COL_PHASE,gap_thr,...
        FS_TARGET,WRAP_ON,DWT_WAVE,DWT_LEVEL,KEEPD,...
        winN,hopN,FMIN_HZ,FMAX_HZ,NFFT_BIG,CONF_THR,FINAL_MOVMEAN_SEC);
    hNC_p=hNC;
    if okNC&&FILL_GAPS_FOR_PLOT,hNC_p=fill_gaps_limited(tP,hNC,MAX_GAP_SEC_PLOT);end

    [R2,M2,P2,C2,n2,~,~,~,~]=calc_metrics_aligned(hNC,hP,tP,30,MAX_LAG_SEC);
    fprintf('%s | RMSE=%.3f MAE=%.3f MAPE=%.2f%% CORR=%.3f N=%d\n',rl,R2,M2,P2,C2,n2);

    ax(d)=nexttile;hold(ax(d),'on');grid(ax(d),'on');
    h1=plot(ax(d),tP,hP,'r-','LineWidth',2);
    h2=plot(ax(d),tP,hNC_p,'k-','LineWidth',1.5);
    ylabel(ax(d),'HR (bpm)','FontSize',11,'FontWeight','bold');
    title(ax(d),ds.name,'FontSize',12,'FontWeight','bold');
    set(ax(d),'FontSize',10);
    if d<numel(DATASETS),set(ax(d),'XTickLabel',[]);else
        xlabel(ax(d),'Tempo (s)','FontSize',11,'FontWeight','bold');end
    if isempty(hPL),hPL=h1;hNCL=h2;end
end
linkaxes(ax,'xy');xlim(ax(1),[0 250]);ylim(ax(1),[60 110]);
lgd=legend([hPL hNCL],{'POLAR H10','HEALTH-RADAR'},...
    'Orientation','horizontal','FontSize',11,'Box','off');
lgd.Layout.Tile='north';
fprintf('\n====== Fim. ======\n');

% =========================================================================
% FUNCOES
% =========================================================================

% Confianca: 3 fatores (sem jump)
function conf=calc_confidence(Pb)
    [pmax,~]=max(Pb);
    f1 = pmax/(sum(Pb)+eps);                           % pico/soma
    f2 = min(pmax/(median(Pb)+eps)/10.0, 1.0);         % pico/mediana
    Ps = sort(Pb,'descend');
    if numel(Ps)>=2&&Ps(2)>0
        f3 = min((Ps(1)/Ps(2))/3.0, 1.0);             % pico/2o_pico
    else
        f3 = 1.0;
    end
    conf = (f1*f2*f3)^(1/3);
end

% FFT peak HR para uma janela, com confianca
function [hr_bpm,conf]=fft_peak_hr(seg,Fs,winN,fmin,fmax,nfft)
    seg=seg(:)-mean(seg,'omitnan');
    seg=seg.*hann(winN);
    nfft=max(nfft,2^nextpow2(winN));
    X=fft(seg,nfft);P=abs(X).^2;
    nh=floor(nfft/2)+1;P=P(1:nh);
    f=(0:nfft-1)'*(Fs/nfft);f=f(1:nh);
    mask=(f>=fmin)&(f<=min(fmax,Fs/2));
    Pb=P(mask);fb=f(mask);
    hr_bpm=NaN;conf=0;
    if numel(Pb)<3||max(Pb)<=0,return;end
    conf=calc_confidence(Pb);
    [~,kmax]=max(Pb);
    if kmax>1&&kmax<numel(Pb)
        y1=log(max(Pb(kmax-1),eps));y2=log(max(Pb(kmax),eps));y3=log(max(Pb(kmax+1),eps));
        den=y1-2*y2+y3;
        if abs(den)>eps,delta=0.5*(y1-y3)/den;delta=max(min(delta,1),-1);
        else,delta=0;end
        fpk=fb(kmax)+delta*(fb(2)-fb(1));
    else,fpk=fb(kmax);end
    hr_bpm=60*fpk;
end

% Non-causal: DWT batch + FFT
function [h_onP,ok]=algo_noncausal(rpath,tP,COL_T,COL_PH,gap_thr,...
    Fs,WRAP,DWT_W,DWT_L,KEEPD,winN,hopN,FMIN,FMAX,NFFT,CTHR,MM_SEC)
    ok=false;h_onP=nan(size(tP));
    try A=readmatrix(rpath);catch,return;end
    if size(A,2)<max(COL_T,COL_PH),return;end
    tp=A(:,COL_T);ph=A(:,COL_PH);tp=tp(:);ph=ph(:);
    m=isfinite(tp)&isfinite(ph)&(abs(ph)<=10);
    tp=tp(m);ph=ph(m);
    if numel(tp)<8,return;end
    t0=(tp-tp(1))/1000;
    [t0,ia]=unique(t0,'stable');ph=ph(ia);
    ph=force_fin(double(ph));ph=unwrap(ph);
    if WRAP==1,x0=wrapph(ph);else,x0=ph-mean(ph,'omitnan');x0=force_fin(x0);end
    segs=seg_gaps(t0,gap_thr);dt=1/Fs;
    ts_l={};hr_l={};
    for s=1:size(segs,1)
        i0=segs(s,1);i1=segs(s,2);
        ts=t0(i0:i1);xs=x0(i0:i1);
        [ts,iu]=unique(ts,'stable');xs=xs(iu);
        if numel(ts)<8,continue;end
        tn=(ts(1):dt:ts(end))';
        xn=interp1(ts,xs,tn,'linear');
        v=isfinite(tn)&isfinite(xn);tn=tn(v);xn=xn(v);
        if numel(tn)<winN,continue;end
        xn=force_fin(xn);
        xd=[0;diff(xn)];xd=force_fin(xd);xd=xd-mean(xd,'omitnan');
        [C,L]=wavedec(xd,DWT_L,DWT_W);
        xr=dwt_keep(C,L,DWT_W,DWT_L,KEEPD);
        N=numel(xr);
        idx0=1:hopN:(N-winN+1);nt=numel(idx0);
        t_c=zeros(nt,1);hr_v=nan(nt,1);
        for i=1:nt
            i1=idx0(i);
            seg=xr(i1:i1+winN-1);
            [hval,cval]=fft_peak_hr(seg,Fs,winN,FMIN,FMAX,NFFT);
            if cval>=CTHR, hr_v(i)=hval; end
            t_c(i)=tn(i1+floor(winN/2));
        end
        if all(isnan(hr_v)),continue;end
        if MM_SEC>0
            W=max(1,round(MM_SEC/(hopN/Fs)));
            hr_v=movmean(hr_v,W,'omitnan','Endpoints','shrink');
        end
        ts_l{end+1}=t_c(:);hr_l{end+1}=hr_v(:);
    end
    [tR,hR]=cat_segs(ts_l,hr_l);
    if isempty(tR),return;end
    v=isfinite(tR)&isfinite(hR);
    if nnz(v)<5,return;end
    h_onP=interp1(tR(v),hR(v),tP,'linear',NaN);ok=true;
end

% Real-time causal: DWT sobre buffer + FFT com NFFT grande
function [t_out,hr_out]=algo_realtime(A_raw,COL_T,COL_PH,gap_thr,...
    Fs,WRAP,DWT_W,DWT_L,KEEPD,winN,hopN,FMIN,FMAX,NFFT,CTHR,MM_SEC)
    t_out=[];hr_out=[];
    tp=A_raw(:,COL_T);ph=A_raw(:,COL_PH);
    m=isfinite(tp)&isfinite(ph)&(abs(ph)<=10);
    tp=tp(m);ph=ph(m);
    if numel(tp)<8,return;end
    t0=(tp-tp(1))/1000;
    [t0,ia]=unique(t0,'stable');ph=ph(ia);
    phu=unwrap(force_fin(double(ph)));
    if WRAP==1,x0=wrapph(phu);else,x0=phu-mean(phu);end
    dt=1/Fs;
    buf_x=zeros(winN,1);buf_t=zeros(winN,1);
    buf_n=0;hop_cnt=0;
    HOP_S=hopN/Fs;W_mm=max(1,round(MM_SEC/HOP_S));
    hr_hist=nan(W_mm,1);
    t_tgt=t0(1);ri=2;
    while ri<=length(t0)
        tp_=t0(ri-1);xp=x0(ri-1);
        tn_=t0(ri);xn=x0(ri);
        if (tn_-tp_)>gap_thr
            buf_n=0;hop_cnt=0;hr_hist=nan(W_mm,1);
            buf_x=zeros(winN,1);buf_t=zeros(winN,1);
            t_tgt=tn_;ri=ri+1;continue;
        end
        while t_tgt<=tn_
            if t_tgt>=tp_
                xi=xp+(xn-xp)*((t_tgt-tp_)/(tn_-tp_));
                if buf_n<winN
                    buf_n=buf_n+1;buf_x(buf_n)=xi;buf_t(buf_n)=t_tgt;
                else
                    buf_x=[buf_x(2:end);xi];buf_t=[buf_t(2:end);t_tgt];
                    hop_cnt=hop_cnt+1;
                end
                if buf_n==winN&&hop_cnt>=hopN
                    hop_cnt=0;
                    xp_=buf_x-mean(buf_x,'omitnan');
                    xd=[0;diff(xp_)];xd=force_fin(xd);xd=xd-mean(xd,'omitnan');
                    [C,L]=wavedec(xd,DWT_L,DWT_W);
                    xhr=dwt_keep(C,L,DWT_W,DWT_L,KEEPD);
                    [hr_val,cval]=fft_peak_hr(xhr,Fs,winN,FMIN,FMAX,NFFT);
                    if cval<CTHR, hr_val=NaN; end
                    if MM_SEC>0
                        hr_hist=[hr_hist(2:end);hr_val];
                        hr_final=mean(hr_hist,'omitnan');
                    else,hr_final=hr_val;end
                    t_out=[t_out;buf_t(end)];
                    hr_out=[hr_out;hr_final];
                end
            end
            t_tgt=t_tgt+dt;
        end
        ri=ri+1;
    end
end

% DWT reconstrucao parcial
function xr=dwt_keep(C,L,wname,Nlev,keepD)
    keepD=unique(keepD(:))';keepD=keepD(keepD>=1&keepD<=Nlev);
    Cn=zeros(size(C));idx=1+L(1);
    for j=1:Nlev
        lev=Nlev-j+1;lenD=L(j+1);i1=idx;i2=idx+lenD-1;
        if any(keepD==lev),Cn(i1:i2)=C(i1:i2);end
        idx=i2+1;
    end
    xr=waverec(Cn,L,wname);xr=xr(:);
    xr=xr-mean(xr,'omitnan');
    xr=fillmissing(xr,'linear','EndValues','nearest');
    xr(~isfinite(xr))=0;
end

% Metricas
function [RMSE,MAE,MAPE,CORR,nV]=calc_metrics(h_est,h_ref,min_n_corr)
    m=isfinite(h_est)&isfinite(h_ref);nV=nnz(m);
    if nV<5,RMSE=NaN;MAE=NaN;MAPE=NaN;CORR=NaN;return;end
    e=h_est(m)-h_ref(m);
    RMSE=sqrt(mean(e.^2,'omitnan'));
    MAE=mean(abs(e),'omitnan');
    MAPE=mean(abs(e./h_ref(m)),'omitnan')*100;
    if nV>=min_n_corr,CORR=corr(h_est(m),h_ref(m),'Rows','complete');
    else,CORR=NaN;end
end

% Metricas com alinhamento por xcorr
function [RMSE,MAE,MAPE,CORR,nV,lag_opt,rmax,lag_sec,dt_eff]=calc_metrics_aligned(h_est,h_ref,tP,min_n,max_lag_s)
    lag_opt=0;rmax=NaN;lag_sec=0;dt_eff=NaN;
    m=isfinite(h_est)&isfinite(h_ref);
    if nnz(m)<min_n
        [RMSE,MAE,MAPE,CORR,nV]=calc_metrics(h_est,h_ref,5);return;end
    [i0,i1]=longest_run(m);
    if (i1-i0+1)<min_n
        [RMSE,MAE,MAPE,CORR,nV]=calc_metrics(h_est,h_ref,5);return;end
    t_run=tP(i0:i1);
    dt_eff=(t_run(end)-t_run(1))/(numel(t_run)-1);
    x=h_est(i0:i1)-mean(h_est(i0:i1));
    y=h_ref(i0:i1)-mean(h_ref(i0:i1));
    [r,lags]=xcorr(x,y,'coeff');
    max_lag_samp=round(max_lag_s/dt_eff);
    r_s=r;r_s(abs(lags)>max_lag_samp)=-Inf;
    [rmax,ix]=max(r_s);
    lag_opt=lags(ix);lag_sec=lag_opt*dt_eff;
    N=numel(h_est);h_al=nan(N,1);
    if lag_opt>=0
        if lag_opt<N,h_al(1:N-lag_opt)=h_est(lag_opt+1:N);end
    else
        al=abs(lag_opt);
        if al<N,h_al(al+1:N)=h_est(1:N-al);end
    end
    [RMSE,MAE,MAPE,CORR,nV]=calc_metrics(h_al,h_ref,5);
end

function [i0,i1]=longest_run(m)
    m=m(:);d=diff([0;m;0]);
    s=find(d==1);e=find(d==-1)-1;
    [~,b]=max(e-s+1);i0=s(b);i1=e(b);
end

% Utilitarias
function segs=seg_gaps(t,gap)
    t=t(:);n=numel(t);
    if n<2,segs=[1 n];return;end
    brk=find(diff(t)>gap);
    if isempty(brk),segs=[1 n];else,segs=[[1;brk+1],[brk;n]];end
end

function [ta,ha]=cat_segs(tl,hl)
    ta=[];ha=[];
    for k=1:numel(tl)
        ts=tl{k}(:);hs=hl{k}(:);
        m=isfinite(ts)&isfinite(hs);ts=ts(m);hs=hs(m);
        if numel(ts)<2,continue;end
        [ts,iu]=unique(ts,'stable');hs=hs(iu);
        ta=[ta;ts];ha=[ha;hs];
    end
    if isempty(ta),return;end
    [ta,iu]=unique(ta,'stable');ha=ha(iu);
end

function x=wrapph(ph)
    ph=double(ph(:));ph=force_fin(ph);
    x=mod(ph+pi,2*pi)-pi;x=x-mean(x,'omitnan');
end

function x=force_fin(x)
    x=double(x(:));ok=isfinite(x);
    if all(ok),return;end
    if nnz(ok)>=2
        xi=find(ok);x(~ok)=interp1(xi,x(ok),find(~ok),'linear','extrap');
    elseif nnz(ok)==1,x(~ok)=x(ok);
    else,x(:)=0;return;end
    x(~isfinite(x))=0;
end

function p2=resolve_any_path(p)
    p=string(p);
    if exist(p,'file')==2,p2=char(p);return;end
    exts=[".csv" ".txt" ".CSV" ".TXT"];
    for e=exts,if exist(p+e,'file')==2,p2=char(p+e);return;end,end
    p2=char(p);
end

function s=get_file_stem(p)
    [~,s,~]=fileparts(char(p));s=string(s);
end

function [t_sec,HR]=read_txt_polar_flex(p)
    fid=fopen(p,'r');
    if fid==-1,t_sec=[];HR=[];return;end
    raw=textscan(fid,'%s','Delimiter','\n');raw=raw{1};fclose(fid);
    if numel(raw)<2,t_sec=[];HR=[];return;end
    lines=raw(2:end);
    ts_str=strings(numel(lines),1);HR=nan(numel(lines),1);
    for i=1:numel(lines)
        L=string(lines{i});
        if contains(L,";"),parts=split(L,";");else,parts=split(L,",");end
        if numel(parts)>=2
            ts_str(i)=strtrim(parts(1));HR(i)=str2double(strtrim(parts(2)));
        end
    end
    try,t_dt=datetime(ts_str,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS");
    catch,t_dt=datetime(ts_str,'InputFormat',"yyyy-MM-dd HH:mm:ss.SSS");end
    t_sec=seconds(t_dt-t_dt(1));t_sec=t_sec(:);HR=HR(:);
    m=isfinite(t_sec)&isfinite(HR);t_sec=t_sec(m);HR=HR(m);
    [t_sec,iu]=unique(t_sec,'stable');HR=HR(iu);
end

function yp=fill_gaps_limited(t,y,max_gap)
    t=t(:);y=y(:);yp=y;
    if numel(t)<3,return;end
    isn=~isfinite(yp);if ~any(isn),return;end
    d=diff([0;isn;0]);iS=find(d==1);iE=find(d==-1)-1;
    for k=1:numel(iS)
        a=iS(k);b=iE(k);
        if a==1||b==numel(yp),continue;end
        if (t(b)-t(a))<=max_gap
            yp(a:b)=interp1([t(a-1);t(b+1)],[yp(a-1);yp(b+1)],t(a:b),'linear');
        end
    end
end