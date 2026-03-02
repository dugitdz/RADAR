clc; clear; close all;

BASE = "C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\";

DATASETS = [
    struct("name","H2",    "radar", BASE+"phases.csv",      "polar", BASE+"POLARH2.txt")
    struct("name","H3",    "radar", BASE+"phases_raw.csv",  "polar", BASE+"POLARH3.txt")
    struct("name","H10",   "radar", BASE+"phase.csv",       "polar", BASE+"POLARH10.txt")
    struct("name","TESTE", "radar", BASE+"Radar_teste.csv", "polar", BASE+"Polar_teste.txt")
];

COL_T_MS=1; COL_PHASE=2;

gap_thr=0.5;
FS=15.0;
IMETH='linear';

DWT_WAVE='db5';
DWT_LEVEL=4;
HR_LEAF="D4D3";          % "D3" ou "D4D3" (seu atalho)

WIN_N=128;             % janela HR (amostras)
HOP_N=32;              % hop HR (amostras)  (20 => 1s em Fs=20)
FMIN_HZ=0.8;
FMAX_HZ=3.0;

CONF_MODE="peak_over_median";
CONF_THR=1.5;

FINAL_MOVMEAN_SEC = 7; % movmean EM SEGUNDOS

SPEC_BPM_MAX=300;      % só pro plot de FFT das folhas (single FFT)

fprintf('\n=== MRVS sem AKF | somente janela RECT ===\n');
fprintf('Fs=%.2f | WIN_N=%d HOP_N=%d | band=[%.2f %.2f] Hz | movmean=%gs | leaf=%s | wave=%s L%d\n', ...
    FS, WIN_N, HOP_N, FMIN_HZ, FMAX_HZ, FINAL_MOVMEAN_SEC, HR_LEAF, DWT_WAVE, DWT_LEVEL);
fprintf('=============================================================\n');

fig = figure('Color','w','Name','MRVS sem AKF | HR vs Polar (RECT)');
tg  = uitabgroup(fig);

for d=1:numel(DATASETS)
    ds = DATASETS(d);

    A = readmatrix(ds.radar);
    tp = A(:,COL_T_MS); ph0 = A(:,COL_PHASE);
    tp = tp(:); ph0 = ph0(:);

    t0 = (tp - tp(1))/1000;
    [t0, ia] = unique(t0,'stable');
    ph0 = ph0(ia);

    ph = force_finite_vector(double(ph0));
    ph = unwrap(ph);

    segs = segment_by_gaps(t0, gap_thr);

    [tP,hP] = read_txt_polar_flex(ds.polar);
    if isempty(tP)
        fprintf('%s | polar vazio\n', ds.name);
        continue;
    end
    tP=tP(:); hP=hP(:);

    fprintf('\n---- %s ----\n', ds.name);

    tab = uitab(tg,'Title',ds.name);
    ax = axes(tab); hold(ax,'on'); grid(ax,'on');
    plot(ax,tP,hP,'c--','LineWidth',1.6);

    % acumula HR por segmentos
    tseg = {}; hrseg = {};

    % maior segmento (xdiff) para plot das folhas
    best_len=-inf; best_xdiff=[];

    dt = 1/FS;

    for s=1:size(segs,1)
        i0=segs(s,1); i1=segs(s,2);
        ts=t0(i0:i1); ps=ph(i0:i1);

        [ts,iu]=unique(ts,'stable'); ps=ps(iu);
        if numel(ts)<8, continue; end

        tnew = (ts(1):dt:ts(end))';
        pnew = interp1(ts, ps, tnew, IMETH);

        ok=isfinite(tnew)&isfinite(pnew);
        tnew=tnew(ok); pnew=pnew(ok);
        if numel(tnew)<(WIN_N+2), continue; end

        pnew = force_finite_vector(pnew);

        % MRVS: diferenciação da fase unwrapped
        xdiff = [0; diff(pnew)];
        xdiff = force_finite_vector(xdiff);
        xdiff = xdiff - mean(xdiff,'omitnan');

        seg_len = tnew(end)-tnew(1);
        if seg_len > best_len
            best_len  = seg_len;
            best_xdiff = xdiff(:);
        end

        % DWT -> leaf (SEM AKF)
        leaf = dwt_leaf_from_xdiff(xdiff, DWT_WAVE, DWT_LEVEL, HR_LEAF);

        % (mantido) normalização por energia (não é AKF)
        leaf_ref = leaf;
        leaf = sqrt_energy_match(leaf, leaf_ref);

        % HR por FFT (RECT)
        [t_ctr, hr_bpm] = estimate_hr_fftN_rect(tnew, leaf, FS, WIN_N, HOP_N, FMIN_HZ, FMAX_HZ, CONF_MODE, CONF_THR);
        if isempty(t_ctr), continue; end

        % movmean em segundos
        if FINAL_MOVMEAN_SEC > 0
            W = max(1, round(FINAL_MOVMEAN_SEC / (HOP_N/FS)));
            hr_bpm = movmean(hr_bpm, W, 'Endpoints','shrink');
        end

        tseg{end+1}=t_ctr(:);
        hrseg{end+1}=hr_bpm(:);
    end

    [tR,hR]=concat_segments(tseg,hrseg);
    if isempty(tR)
        title(ax, sprintf('%s | sem HR válido', ds.name));
        xlabel(ax,'t (s)'); ylabel(ax,'HR (bpm)');
        continue;
    end

    h_onP = interp1(tR,hR,tP,'linear',NaN);

    [RMSE,MAE,CORR,Nv]=metrics_on_mask(h_onP,hP);
    fprintf('rect | RMSE=%7.3f MAE=%7.3f CORR=%7.3f N=%d\n', RMSE, MAE, CORR, Nv);

    plot(ax, tP, h_onP, 'k-', 'LineWidth', 1.3);
    legend(ax, {'POLAR','RADAR(rect)'}, 'Location','best');
    title(ax, sprintf('%s | rect | RMSE=%.3f MAE=%.3f Corr=%.3f | N=%d', ds.name, RMSE, MAE, CORR, Nv));
    xlabel(ax,'t (s)'); ylabel(ax,'HR (bpm)');

    % FFT das folhas (janela separada por dataset)
    if ~isempty(best_xdiff) && numel(best_xdiff) >= 64
        plot_leaves_fft_single_rect(ds.name, best_xdiff, FS, DWT_WAVE, DWT_LEVEL, HR_LEAF, SPEC_BPM_MAX);
    end
end

fprintf('\nDone.\n');

% ===================== FUNÇÕES =====================

function leaf = dwt_leaf_from_xdiff(xdiff, wname, level, HR_LEAF)
    [C,L]=wavedec(xdiff,level,wname);
    D=cell(level,1);
    for lev=1:level
        D{lev}=wrcoef('d',C,L,wname,lev);
    end
    if HR_LEAF=="D4D3"
        leaf = D{4}+D{3};
    else
        leaf = D{3};
    end
    leaf = force_finite_vector(leaf);
end

function [t_ctr, hr_bpm] = estimate_hr_fftN_rect(t, x, Fs, winN, hopN, fmin, fmax, conf_mode, conf_thr)
    t=t(:); x=double(x(:));
    x=x-mean(x,'omitnan');
    x=fillmissing(x,'linear','EndValues','nearest'); x(~isfinite(x))=0;

    N=numel(x);
    if N<winN, t_ctr=[]; hr_bpm=[]; return; end
    hopN=max(1,min(hopN,winN));

    nfft=2^nextpow2(winN);
    nh=floor(nfft/2)+1;
    f=(0:nh-1)'*(Fs/nfft);

    band=(f>=fmin)&(f<=min(fmax,Fs/2));
    if nnz(band)<3, t_ctr=[]; hr_bpm=[]; return; end

    idx0=1:hopN:(N-winN+1);
    nt=numel(idx0);
    t_ctr=zeros(nt,1); hr_bpm=nan(nt,1);

    for i=1:nt
        i1=idx0(i);
        seg=x(i1:i1+winN-1);
        seg=seg-mean(seg,'omitnan');   % rect window => sem multiplicar por w

        X=fft(seg,nfft);
        P=abs(X).^2;
        P=P(1:nh);

        Pb=P(band); fb=f(band);
        t_ctr(i)=t(i1+floor(winN/2));

        if ~any(isfinite(Pb)) || max(Pb)<=0, continue; end
        [pmax,k]=max(Pb);

        ok=false;
        if conf_mode=="peak_over_median"
            med=median(Pb(Pb>0)); if ~isfinite(med)||med<=0, med=eps; end
            ok=(pmax/med)>=conf_thr;
        elseif conf_mode=="pmax_norm"
            ok=(pmax/(sum(Pb)+eps))>=conf_thr;
        else
            ok=true;
        end

        if ok, hr_bpm(i)=60*fb(k); end
    end
end

function plot_leaves_fft_single_rect(name, xdiff, Fs, wname, level, HR_LEAF, BPM_MAX)
    [C,L]=wavedec(xdiff,level,wname);
    A = wrcoef('a',C,L,wname,level);
    D = cell(level,1);
    for lev=1:level
        D{lev}=wrcoef('d',C,L,wname,lev);
    end
    if HR_LEAF=="D4D3"
        leaf = D{4}+D{3};
    else
        leaf = D{3};
    end

    nrows = level + 3; % xdiff + A + D(level..1) + leaf
    f=figure('Color','w','Name',sprintf('%s - FFT folhas (RECT) | %s L%d',name,wname,level));
    tl=tiledlayout(f,nrows,1,'Padding','compact','TileSpacing','compact');

    ax=nexttile(tl,1); hold(ax,'on'); grid(ax,'on');
    [bpm,Pdb]=fft_db_single_rect(xdiff,Fs);
    plot(ax,bpm,Pdb,'LineWidth',1.0); xlim(ax,[0 BPM_MAX]); ylim(ax,[-70 5]);
    ylabel(ax,'FFT xdiff'); title(ax,sprintf('%s | leaf=%s',name,HR_LEAF));

    ax=nexttile(tl,2); hold(ax,'on'); grid(ax,'on');
    [bpm,Pdb]=fft_db_single_rect(A,Fs);
    plot(ax,bpm,Pdb,'LineWidth',1.0); xlim(ax,[0 BPM_MAX]); ylim(ax,[-70 5]);
    ylabel(ax,sprintf('FFT A%d',level));

    rr=3;
    for lev=level:-1:1
        ax=nexttile(tl,rr); hold(ax,'on'); grid(ax,'on');
        [bpm,Pdb]=fft_db_single_rect(D{lev},Fs);
        plot(ax,bpm,Pdb,'LineWidth',1.0); xlim(ax,[0 BPM_MAX]); ylim(ax,[-70 5]);
        ylabel(ax,sprintf('FFT D%d',lev));
        rr=rr+1;
    end

    ax=nexttile(tl,nrows); hold(ax,'on'); grid(ax,'on');
    [bpm,Pdb]=fft_db_single_rect(leaf,Fs);
    plot(ax,bpm,Pdb,'LineWidth',1.2); xlim(ax,[0 BPM_MAX]); ylim(ax,[-70 5]);
    ylabel(ax,'FFT LEAF'); xlabel(ax,'BPM');
end

function [bpm,Pdb] = fft_db_single_rect(x, Fs)
    x=double(x(:));
    x=x-mean(x,'omitnan');
    x=fillmissing(x,'linear','EndValues','nearest'); x(~isfinite(x))=0;

    N=numel(x);
    nfft=2^nextpow2(N);
    X=fft(x,nfft);
    P=abs(X).^2;

    f=(0:nfft-1)'*(Fs/nfft);
    nh=floor(nfft/2)+1;
    f=f(1:nh); P=P(1:nh);

    bpm=60*f;
    Pref=max(P)+eps;
    Pdb=10*log10(P/Pref+eps);
end

function y = sqrt_energy_match(x_rec, x_ref)
    x_rec=double(x_rec(:)); x_ref=double(x_ref(:));
    Er=sum(x_rec.^2,'omitnan'); E0=sum(x_ref.^2,'omitnan');
    if ~isfinite(Er)||Er<=0||~isfinite(E0)||E0<=0, y=x_rec; return; end
    y = sqrt(E0/Er) * x_rec;
end

function segments = segment_by_gaps(t, gap_thr)
    t=t(:); n=numel(t);
    if n<2, segments=[1 n]; return; end
    brk=find(diff(t)>gap_thr);
    if isempty(brk), segments=[1 n];
    else, segments=[[1; brk+1],[brk; n]];
    end
end

function [t_all,h_all] = concat_segments(tseg_list,hseg_list)
    t_all=[]; h_all=[];
    for k=1:numel(tseg_list)
        tS=tseg_list{k}; hS=hseg_list{k};
        if numel(tS)<2, continue; end
        [tS,iu]=unique(tS,'stable'); hS=hS(iu);
        t_all=[t_all; tS]; h_all=[h_all; hS];
    end
    [t_all,iu]=unique(t_all,'stable'); h_all=h_all(iu);
end

function x = force_finite_vector(x)
    x=double(x(:));
    ok=isfinite(x);
    if all(ok), return; end
    if nnz(ok)>=2
        xi=find(ok);
        x(~ok)=interp1(xi,x(ok),find(~ok),'linear','extrap');
    elseif nnz(ok)==1
        x(~ok)=x(ok);
    else
        x(:)=0; return;
    end
    x(~isfinite(x))=0;
    if ~any(isfinite(x)), x(:)=0; end
end

function [RMSE,MAE,CORR,nValid] = metrics_on_mask(h_est,h_ref)
    m=isfinite(h_est)&isfinite(h_ref);
    nValid=nnz(m);
    if nValid>=5
        e=h_est(m)-h_ref(m);
        RMSE=sqrt(mean(e.^2,'omitnan'));
        MAE=mean(abs(e),'omitnan');
        CORR=corr(h_est(m),h_ref(m),'Rows','complete');
    else
        RMSE=NaN; MAE=NaN; CORR=NaN;
    end
end

function [t_sec,HR] = read_txt_polar_flex(p)
    fid=fopen(p,'r');
    if fid==-1, t_sec=[]; HR=[]; return; end
    raw=textscan(fid,'%s','Delimiter','\n'); raw=raw{1}; fclose(fid);
    if numel(raw)<2, t_sec=[]; HR=[]; return; end
    lines=raw(2:end);
    ts_str=strings(numel(lines),1);
    HR=nan(numel(lines),1);
    for i=1:numel(lines)
        L=string(lines{i});
        if contains(L,";"), parts=split(L,";"); else, parts=split(L,","); end
        if numel(parts)>=2
            ts_str(i)=strtrim(parts(1));
            HR(i)=str2double(strtrim(parts(2)));
        end
    end
    try
        t_dt=datetime(ts_str,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS");
    catch
        t_dt=datetime(ts_str,'InputFormat',"yyyy-MM-dd HH:mm:ss.SSS");
    end
    t_sec=seconds(t_dt-t_dt(1)); t_sec=t_sec(:); HR=HR(:);
    m=isfinite(t_sec)&isfinite(HR);
    t_sec=t_sec(m); HR=HR(m);
    [t_sec,iu]=unique(t_sec,'stable'); HR=HR(iu);
end