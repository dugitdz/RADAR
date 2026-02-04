// radar_polar_compare.c
// gcc -O2 radar_polar_compare.c -lfftw3 -lm -o radar_polar_compare

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#ifdef _WIN32
#define timegm _mkgmtime
#endif

// ===================== PARÂMETROS CONFIGURÁVEIS ===================== //
static const char *PHASES_CSV_PATH = "C:\\Users\\eduar\\UTFPR\\IC\\PROJ_ESP32\\output\\phases.csv";
static const char *POLAR_TXT_PATH  = "C:\\Users\\eduar\\UTFPR\\IC\\PROJ_ESP32\\output\\POLARH2.txt";

static const int    N        = 16;
static const double overlap  = 0.8;
static const int    win_smooth_raw   = 64;
static const int    win_smooth_final = 64;

static const double gap_thr = 0.5;

static const double dt_grid = 0.1;
static const double IGNORE_START_SEC = 0.0;
// =================================================================== //

static void die(const char *msg){ fprintf(stderr,"ERROR: %s\n",msg); exit(1); }
static int is_nan(double x){ return x != x; }

typedef struct {
    double *tpuro_ms; // col1
    double *breath;   // col3
    double *heart;    // col4
    int n;
} PhaseData;

static int split_line(char *s, char **out, int max_tokens){
    int k=0; char *p=s;
    while (*p && k<max_tokens){
        while (*p==' '||*p=='\t'||*p=='\r'||*p=='\n'||*p==','||*p==';') p++;
        if (!*p) break;
        out[k++]=p;
        while (*p && *p!=',' && *p!=';' && *p!='\t' && *p!='\r' && *p!='\n') p++;
        if (*p){ *p='\0'; p++; }
    }
    return k;
}

static PhaseData read_phases_raw_csv(const char *path){
    FILE *f=fopen(path,"rb");
    if(!f) die("cannot open phases_raw.csv");

    int lines=0;
    char buf[4096];
    while(fgets(buf,sizeof(buf),f)){
        char *p=buf; while(*p==' '||*p=='\t') p++;
        if(*p=='\0'||*p=='\r'||*p=='\n') continue;
        lines++;
    }
    rewind(f);
    if(lines<=0) die("phases_raw.csv empty");

    double *tp=malloc((size_t)lines*sizeof(double));
    double *br=malloc((size_t)lines*sizeof(double));
    double *hr=malloc((size_t)lines*sizeof(double));
    if(!tp||!br||!hr) die("malloc phases arrays");

    int i=0;
    while(fgets(buf,sizeof(buf),f)){
        char *p=buf; while(*p==' '||*p=='\t') p++;
        if(*p=='\0'||*p=='\r'||*p=='\n') continue;

        char *tok[16]={0};
        int ntok=split_line(p,tok,16);
        if(ntok<4) continue;

        tp[i]=strtod(tok[0],NULL);
        br[i]=strtod(tok[2],NULL);
        hr[i]=strtod(tok[3],NULL);
        i++;
    }
    fclose(f);

    if(i<4) die("phases_raw.csv: too few rows");

    PhaseData out={tp,br,hr,i};
    return out;
}

typedef struct { int start,end; } Segment;

static Segment* segment_by_gaps(const double *t,int n,double gap,int *outNSeg,double *out_dt,double *out_srate){
    double *dt_all = (double*)malloc((size_t)(n-1)*sizeof(double));
    if(!dt_all) die("malloc dt_all");

    for(int i=0;i<n-1;i++) dt_all[i]=t[i+1]-t[i];

    int nb=0;
    for(int i=0;i<n-1;i++) if(dt_all[i]>gap) nb++;

    int nSeg = (nb==0)? 1 : (nb+1);
    Segment *segs = (Segment*)malloc((size_t)nSeg*sizeof(Segment));
    if(!segs) die("malloc segs");

    if(nb==0){
        segs[0].start=0; segs[0].end=n-1;
    } else {
        int s=0;
        segs[s].start=0;
        for(int i=0;i<n-1;i++){
            if(dt_all[i]>gap){
                segs[s].end=i;
                s++;
                segs[s].start=i+1;
            }
        }
        segs[s].end=n-1;
    }

    double sum_good=0; int cnt_good=0;
    double sum_pos=0;  int cnt_pos=0;
    for(int i=0;i<n-1;i++){
        if(dt_all[i]>0){ sum_pos += dt_all[i]; cnt_pos++; }
        if(dt_all[i]>0 && dt_all[i]<=gap){ sum_good += dt_all[i]; cnt_good++; }
    }

    double dt;
    if(cnt_good>0) dt = sum_good/(double)cnt_good;
    else if(cnt_pos>0) dt = sum_pos/(double)cnt_pos;
    else die("all dt non-positive");

    *out_dt = dt;
    *out_srate = 1.0/dt;
    *outNSeg = nSeg;

    free(dt_all);
    return segs;
}

typedef struct {
    int L;
    double *in;
    fftw_complex *out;
    fftw_plan plan;
    int ready;
} PlanSlot;

static void ensure_plan(PlanSlot *slots, int L){
    if(L<4) die("ensure_plan L<4");
    if(slots[L].ready) return;
    slots[L].L=L;
    slots[L].in  = (double*)fftw_malloc((size_t)L*sizeof(double));
    slots[L].out = (fftw_complex*)fftw_malloc((size_t)(L/2+1)*sizeof(fftw_complex));
    if(!slots[L].in || !slots[L].out) die("fftw_malloc window");
    slots[L].plan = fftw_plan_dft_r2c_1d(L, slots[L].in, slots[L].out, FFTW_ESTIMATE);
    if(!slots[L].plan) die("fftw_plan window");
    slots[L].ready=1;
}

static void destroy_plans(PlanSlot *slots, int maxL){
    for(int L=0;L<=maxL;L++){
        if(slots[L].ready){
            fftw_destroy_plan(slots[L].plan);
            fftw_free(slots[L].in);
            fftw_free(slots[L].out);
            slots[L].ready=0;
        }
    }
}

static int estimate_bpm_window(PlanSlot *slots,int maxL,const double *tseg,const double *xseg,int L,double *out_bpm){
    if(L<4 || L>maxL) return 0;

    double dt_sum=0;
    for(int i=0;i<L-1;i++){
        double d = tseg[i+1]-tseg[i];
        if(d<=0) return 0;
        dt_sum += d;
    }
    double fs = 1.0/(dt_sum/(double)(L-1));

    double mean=0;
    for(int i=0;i<L;i++) mean += xseg[i];
    mean /= (double)L;

    ensure_plan(slots,L);
    for(int i=0;i<L;i++) slots[L].in[i] = xseg[i]-mean;

    fftw_execute(slots[L].plan);

    int half = L/2;
    if(half<3) return 0;

    int imax=1; double best=0;
    for(int k=1;k<half;k++){ // ignora DC
        double re=slots[L].out[k][0], im=slots[L].out[k][1];
        double mag = sqrt(re*re+im*im);
        if(mag>best){ best=mag; imax=k; }
    }

    double f_est = (double)imax * fs / (double)L;
    *out_bpm = f_est*60.0;
    return 1;
}

static void movmean_omitnan_segment(const double *x,int a,int b,int win,double *out){
    int n=b-a+1; if(n<=0) return;
    int left=(win-1)/2;
    int right=(win-1)-left;

    double *ps=(double*)calloc((size_t)(n+1),sizeof(double));
    int *pc=(int*)calloc((size_t)(n+1),sizeof(int));
    if(!ps||!pc) die("movmean calloc");

    for(int i=0;i<n;i++){
        double v=x[a+i];
        ps[i+1]=ps[i]; pc[i+1]=pc[i];
        if(!is_nan(v)){ ps[i+1]+=v; pc[i+1]++; }
    }

    for(int i=0;i<n;i++){
        int L=i-left; if(L<0) L=0;
        int R=i+right; if(R>n-1) R=n-1;
        double sum = ps[R+1]-ps[L];
        int cnt = pc[R+1]-pc[L];
        out[a+i] = (cnt>0) ? (sum/(double)cnt) : NAN;
    }

    free(ps); free(pc);
}

static void interp_fill_nan_linear_extrap(const double *t,double *y,int n){
    int first=-1,last=-1,cnt=0;
    for(int i=0;i<n;i++){
        if(!is_nan(y[i])){ if(first<0) first=i; last=i; cnt++; }
    }
    if(cnt==0) return;
    if(cnt==1){
        for(int i=0;i<n;i++) y[i]=y[first];
        return;
    }

    int *g=(int*)malloc((size_t)cnt*sizeof(int));
    if(!g) die("malloc g");
    int k=0;
    for(int i=0;i<n;i++) if(!is_nan(y[i])) g[k++]=i;

    { // extrap left
        int i0=g[0], i1=g[1];
        double m=(y[i1]-y[i0])/(t[i1]-t[i0]);
        for(int i=0;i<i0;i++) y[i]=y[i0]+m*(t[i]-t[i0]);
    }
    { // extrap right
        int i0=g[cnt-2], i1=g[cnt-1];
        double m=(y[i1]-y[i0])/(t[i1]-t[i0]);
        for(int i=i1+1;i<n;i++) y[i]=y[i1]+m*(t[i]-t[i1]);
    }
    for(int j=0;j<cnt-1;j++){
        int i0=g[j], i1=g[j+1];
        double m=(y[i1]-y[i0])/(t[i1]-t[i0]);
        for(int i=i0+1;i<i1;i++){
            if(is_nan(y[i])) y[i]=y[i0]+m*(t[i]-t[i0]);
        }
    }

    free(g);
}

typedef struct {
    double *raw;
    double *final;
} FreqOut;

static FreqOut estimate_rate(const double *t,const double *sig,int nt,const Segment *segs,int nSeg){
    int hop = (int)llround((double)N*(1.0-overlap));
    if(hop<1) hop=1;
    int halfN = N/2;

    double *wsum=(double*)calloc((size_t)nt,sizeof(double));
    double *wcnt=(double*)calloc((size_t)nt,sizeof(double));
    double *raw=(double*)malloc((size_t)nt*sizeof(double));
    double *sm1=(double*)malloc((size_t)nt*sizeof(double));
    double *fin=(double*)malloc((size_t)nt*sizeof(double));
    if(!wsum||!wcnt||!raw||!sm1||!fin) die("malloc rate arrays");

    for(int i=0;i<nt;i++){ raw[i]=NAN; sm1[i]=NAN; fin[i]=NAN; }

    int maxL=N;
    PlanSlot *slots=(PlanSlot*)calloc((size_t)(maxL+1),sizeof(PlanSlot));
    if(!slots) die("calloc slots");

    int rel_len=2*halfN;
    int *rel=(int*)malloc((size_t)rel_len*sizeof(int));
    if(!rel) die("malloc rel");
    for(int i=0;i<rel_len;i++) rel[i]=i-halfN;

    double tseg[64], xseg[64];
    int idx[64];

    for(int s=0;s<nSeg;s++){
        int a=segs[s].start, b=segs[s].end;
        for(int center=a; center<=b; center+=hop){
            int L=0;
            for(int r=0;r<rel_len;r++){
                int ii=center+rel[r];
                if(ii>=a && ii<=b){
                    idx[L++]=ii;
                    if(L>=64) break;
                }
            }
            if(L<4) continue;

            for(int j=0;j<L;j++){
                int ii=idx[j];
                tseg[j]=t[ii];
                xseg[j]=sig[ii];
            }

            double bpm;
            if(!estimate_bpm_window(slots,maxL,tseg,xseg,L,&bpm)) continue;

            for(int j=0;j<L;j++){
                int ii=idx[j];
                int dist=abs(ii-center);
                double w=(double)((halfN+1)-dist);
                if(w<0) w=0;
                wsum[ii]+=bpm*w;
                wcnt[ii]+=w;
            }
        }
    }

    for(int i=0;i<nt;i++) if(wcnt[i]>0) raw[i]=wsum[i]/wcnt[i];

    for(int s=0;s<nSeg;s++)
        movmean_omitnan_segment(raw,segs[s].start,segs[s].end,win_smooth_raw,sm1);

    interp_fill_nan_linear_extrap(t,sm1,nt);

    for(int s=0;s<nSeg;s++)
        movmean_omitnan_segment(sm1,segs[s].start,segs[s].end,win_smooth_final,fin);

    free(wsum); free(wcnt); free(sm1);
    free(rel);
    destroy_plans(slots,maxL);
    free(slots);

    FreqOut out={raw,fin};
    return out;
}

static int parse_polar_datetime(const char *s, double *out_sec_abs){
    int Y,M,D,h,m,sec,ms;
    if(sscanf(s,"%d-%d-%dT%d:%d:%d.%d",&Y,&M,&D,&h,&m,&sec,&ms)!=7){
        if(sscanf(s,"%d-%d-%d %d:%d:%d.%d",&Y,&M,&D,&h,&m,&sec,&ms)!=7) return 0;
    }
    struct tm tmv; memset(&tmv,0,sizeof(tmv));
    tmv.tm_year=Y-1900; tmv.tm_mon=M-1; tmv.tm_mday=D;
    tmv.tm_hour=h; tmv.tm_min=m; tmv.tm_sec=sec;
    time_t tt=timegm(&tmv);
    if(tt==(time_t)-1) return 0;
    *out_sec_abs=(double)tt + (double)ms/1000.0;
    return 1;
}

typedef struct {
    double *t_sec;
    double *hr;
    int n;
} PolarData;

static PolarData read_polar_txt(const char *path){
    FILE *f=fopen(path,"rb");
    if(!f) die("cannot open polar txt");

    char buf[4096];
    if(!fgets(buf,sizeof(buf),f)) die("polar empty");

    int cap=4096,n=0;
    double *t=(double*)malloc((size_t)cap*sizeof(double));
    double *hr=(double*)malloc((size_t)cap*sizeof(double));
    if(!t||!hr) die("malloc polar arrays");

    double t0_abs=0; int got0=0;

    while(fgets(buf,sizeof(buf),f)){
        char *p=buf; while(*p==' '||*p=='\t') p++;
        if(*p=='\0'||*p=='\r'||*p=='\n') continue;

        char *ts=strtok(p,";");
        char *hrtok=strtok(NULL,";");
        if(!ts||!hrtok) continue;

        while(*ts==' '||*ts=='\t') ts++;
        char *e=ts+strlen(ts);
        while(e>ts && (e[-1]==' '||e[-1]=='\t'||e[-1]=='\r'||e[-1]=='\n')) e--;
        *e='\0';

        double t_abs;
        if(!parse_polar_datetime(ts,&t_abs)) continue;

        double vhr=strtod(hrtok,NULL);
        if(!got0){ t0_abs=t_abs; got0=1; }
        double trel=t_abs-t0_abs;

        if(n>=cap){
            cap*=2;
            t=(double*)realloc(t,(size_t)cap*sizeof(double));
            hr=(double*)realloc(hr,(size_t)cap*sizeof(double));
            if(!t||!hr) die("realloc polar arrays");
        }
        t[n]=trel; hr[n]=vhr; n++;
    }
    fclose(f);
    if(n<2) die("polar too few rows");

    // unique stable (dup exata de tempo)
    double *t2=(double*)malloc((size_t)n*sizeof(double));
    double *h2=(double*)malloc((size_t)n*sizeof(double));
    if(!t2||!h2) die("malloc unique");
    int m=0;
    for(int i=0;i<n;i++){
        if(m==0 || t[i]!=t2[m-1]){ t2[m]=t[i]; h2[m]=hr[i]; m++; }
    }
    free(t); free(hr);

    PolarData out={t2,h2,m};
    return out;
}

static double interp1_linear_extrap_one(const double *tx,const double *yx,int n,double tq){
    if(n<2) return NAN;
    if(tq<=tx[0]){
        double m=(yx[1]-yx[0])/(tx[1]-tx[0]);
        return yx[0]+m*(tq-tx[0]);
    }
    if(tq>=tx[n-1]){
        double m=(yx[n-1]-yx[n-2])/(tx[n-1]-tx[n-2]);
        return yx[n-1]+m*(tq-tx[n-1]);
    }
    int i=0;
    while(i+1<n && !(tx[i]<=tq && tq<=tx[i+1])) i++;
    if(i+1>=n) return NAN;
    double a=(tq-tx[i])/(tx[i+1]-tx[i]);
    return yx[i]+a*(yx[i+1]-yx[i]);
}

static double pearson_corr(const double *x,const double *y,int n){
    double sx=0,sy=0; int c=0;
    for(int i=0;i<n;i++) if(!is_nan(x[i]) && !is_nan(y[i])){ sx+=x[i]; sy+=y[i]; c++; }
    if(c<2) return NAN;
    double mx=sx/(double)c, my=sy/(double)c;
    double sxx=0,syy=0,sxy=0;
    for(int i=0;i<n;i++) if(!is_nan(x[i]) && !is_nan(y[i])){
        double dx=x[i]-mx, dy=y[i]-my;
        sxx+=dx*dx; syy+=dy*dy; sxy+=dx*dy;
    }
    if(sxx<=0 || syy<=0) return NAN;
    return sxy/sqrt(sxx*syy);
}

int main(int argc,char **argv){
    const char *phases_path = PHASES_CSV_PATH;
    const char *polar_path  = POLAR_TXT_PATH;
    if(argc>=3){ phases_path=argv[1]; polar_path=argv[2]; }

    PhaseData pd=read_phases_raw_csv(phases_path);
    int nt=pd.n;

    double *t=(double*)malloc((size_t)nt*sizeof(double));
    if(!t) die("malloc t");
    double t0=pd.tpuro_ms[0];
    for(int i=0;i<nt;i++) t[i]=(pd.tpuro_ms[i]-t0)/1000.0;

    int nSeg=0; double dt=0,srate=0;
    Segment *segs=segment_by_gaps(t,nt,gap_thr,&nSeg,&dt,&srate);

    FreqOut HR = estimate_rate(t,pd.heart,nt,segs,nSeg);
    // Se quiser RR também (sem comparação com Polar), descomenta:
    // FreqOut RR = estimate_rate(t,pd.breath,nt,segs,nSeg);

    PolarData pol = read_polar_txt(polar_path);

    double sum=0; int cnt=0;
    for(int i=0;i<pol.n;i++) if(!is_nan(pol.hr[i])){ sum+=pol.hr[i]; cnt++; }
    double HR_polar_mean = (cnt>0)? (sum/(double)cnt) : NAN;

    double t_start = (t[0] > pol.t_sec[0]) ? t[0] : pol.t_sec[0];
    double t_end   = (t[nt-1] < pol.t_sec[pol.n-1]) ? t[nt-1] : pol.t_sec[pol.n-1];
    if(t_end <= t_start) die("no overlap time between radar and polar");

    int ngrid = (int)floor((t_end - t_start)/dt_grid) + 1;
    double *ref=(double*)malloc((size_t)ngrid*sizeof(double));
    double *test=(double*)malloc((size_t)ngrid*sizeof(double));
    double *tt=(double*)malloc((size_t)ngrid*sizeof(double));
    if(!ref||!test||!tt) die("malloc grid arrays");

    for(int i=0;i<ngrid;i++){
        double tq = t_start + dt_grid*(double)i;
        tt[i]=tq;
        ref[i]  = interp1_linear_extrap_one(pol.t_sec, pol.hr, pol.n, tq);
        test[i] = interp1_linear_extrap_one(t, HR.final, nt, tq);
    }

    int n2=0;
    for(int i=0;i<ngrid;i++) if(tt[i] >= t_start + IGNORE_START_SEC) n2++;
    if(n2<2) die("not enough samples after IGNORE_START_SEC");

    double mae_sum=0, mse_sum=0; int c2=0;
    for(int i=0;i<ngrid;i++){
        if(tt[i] < t_start + IGNORE_START_SEC) continue;
        if(is_nan(ref[i]) || is_nan(test[i])) continue;
        double e = test[i]-ref[i];
        mae_sum += fabs(e);
        mse_sum += e*e;
        c2++;
    }
    if(c2<2) die("not enough valid points for metrics");

    double MAE  = mae_sum/(double)c2;
    double RMSE = sqrt(mse_sum/(double)c2);
    double R = pearson_corr(test, ref, ngrid);

    printf("dt=%.6f s | srate=%.3f Hz | nSeg=%d | hop=%d\n", dt, srate, nSeg, (int)llround((double)N*(1.0-overlap)));
    printf("HR Polar medio: %.4f bpm\n", HR_polar_mean);
    printf("MAE : %.4f bpm\n", MAE);
    printf("RMSE: %.4f bpm\n", RMSE);
    printf("Corr: %.4f\n", R);

    free(tt); free(ref); free(test);
    free(pol.t_sec); free(pol.hr);
    free(HR.raw); free(HR.final);
    // free(RR.raw); free(RR.final);

    free(segs);
    free(t);
    free(pd.tpuro_ms); free(pd.breath); free(pd.heart);
    return 0;
}
