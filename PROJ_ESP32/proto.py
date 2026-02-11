import numpy as np
import pandas as pd
from collections import deque
from scipy.signal import cheby2, sosfilt, sosfilt_zi

# ========================== CONFIG ==========================
BASE = r"C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\\"
DATASETS = [
    {"name": "H2",  "radar": BASE + "phases.csv",     "polar": BASE + "POLARH2.txt"},
    {"name": "H3",  "radar": BASE + "phases_raw.csv", "polar": BASE + "POLARH3.txt"},
    {"name": "H10", "radar": BASE + "phase.csv",      "polar": BASE + "POLARH10.txt"},
]

FS_TARGET, GAP_THR, DT_GRID, F_MIN_HZ = 50.0, 0.6, 0.1, 0.05

CFG_REALTIME = {
    "N": 512, "hop": 64,
    "FILT_ORD": 4, "FILT_WI": 0.595, "FILT_WF": 3.05, "FILT_RS": 30.50,
    "PRIOR_INIT_BPM": 80.0, "PRIOR_MIN_BPM": 35.0, "PRIOR_MAX_BPM": 240.0,
    "PRIOR_SIGMA_HZ": 0.40, "PRIOR_ALPHA": 0.05,
    "RT_EMA_TAU_SEC": 4.0,
}

# ========================== IO & PREP ==========================
def read_radar_pandas(path):
    try:
        df = pd.read_csv(path, sep=r'[;,\t| ]+', engine='python', header=None, comment='#')
        df = df.apply(pd.to_numeric, errors='coerce').dropna()
        if df.shape[1] < 2: return np.array([]), np.array([])
        return df.iloc[:, 0].values.astype(float), df.iloc[:, 1].values.astype(float)
    except Exception:
        return np.array([]), np.array([])

def read_polar_pandas(path):
    try:
        df = pd.read_csv(path, sep=";", header=0, usecols=[0, 1], engine="python")
        df.columns = ['ts', 'hr']
        df['hr'] = pd.to_numeric(df['hr'], errors='coerce')
        df.dropna(subset=['hr', 'ts'], inplace=True)
        if df.empty: return np.array([]), np.array([])
        
        # Pandas detecta formatos de data automaticamente muito melhor que loops manuais
        t_dt = pd.to_datetime(df['ts'].astype(str).str.replace("Z", ""), errors='coerce')
        valid = t_dt.notna()
        t_sec = (t_dt[valid] - t_dt[valid].iloc[0]).dt.total_seconds().values
        hr = df['hr'][valid].values
        
        # Remove duplicatas de tempo e ordena
        _, idx = np.unique(t_sec, return_index=True)
        return t_sec[idx], hr[idx]
    except Exception:
        return np.array([]), np.array([])

def segment_and_resample(t, x, fs_target, gap_thr):
    t, x = np.atleast_1d(t), np.atleast_1d(x)
    if t.size < 2: return np.array([]), np.array([])
    
    # Detecção de gaps
    dt_seq = np.diff(t)
    breaks = np.flatnonzero(dt_seq > gap_thr)
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks, t.size - 1]
    
    # Resample
    dt = 1.0 / fs_target
    t_out, x_out = [], []
    
    for i0, i1 in zip(starts, ends):
        ts, xs = t[i0:i1+1], x[i0:i1+1]
        _, u_idx = np.unique(ts, return_index=True) # remove duplicados locais
        ts, xs = ts[u_idx], xs[u_idx]
        
        if ts.size < 4: continue
        
        tn = np.arange(ts[0], ts[-1] + 1e-12, dt)
        xn = np.interp(tn, ts, xs)
        t_out.append(tn)
        x_out.append(xn)
        
    if not t_out: return np.array([]), np.array([])
    return np.concatenate(t_out), np.concatenate(x_out)

# ========================== METRICS & MATH ==========================
def design_sos_cheby2(fs, order, wi, wf, rs):
    if wf >= fs/2 or wi <= 0 or wi >= wf or (order % 2) != 0: return None, None
    sos = cheby2(order // 2, rs, [wi/(fs/2), wf/(fs/2)], btype="bandpass", output="sos")
    return sos, sosfilt_zi(sos)

def parabolic_peak(P, k):
    if k <= 0 or k >= (len(P) - 1): return 0.0
    val = np.log(P[k-1:k+2] + 1e-20)
    den = val[0] - 2*val[1] + val[2]
    if den == 0 or not np.isfinite(den): return 0.0
    return float(np.clip(0.5 * (val[0] - val[2]) / den, -0.75, 0.75))

def calc_metrics(yhat, yref, t_hat, t_ref):
    # Interpola referência para o grid da estimativa (com validação temporal)
    t_min, t_max = max(t_hat[0], t_ref[0]), min(t_hat[-1], t_ref[-1])
    if t_max <= t_min: return np.nan, np.nan, np.nan, 0
    
    t_common = np.arange(t_min, t_max + 1e-12, DT_GRID)
    yh = np.interp(t_common, t_hat, yhat)
    
    # Garante unicidade da referência antes de interpolar
    _, u_ref = np.unique(t_ref, return_index=True)
    yr = np.interp(t_common, t_ref[u_ref], yref[u_ref])
    
    mask = np.isfinite(yh) & np.isfinite(yr)
    if (n := mask.sum()) < 3: return np.nan, np.nan, np.nan, int(n)
    
    e = yh[mask] - yr[mask]
    mae, rmse = np.mean(np.abs(e)), np.sqrt(np.mean(e**2))
    r = np.corrcoef(yh[mask], yr[mask])[0, 1] if np.std(yh[mask]) > 0 and np.std(yr[mask]) > 0 else np.nan
    return mae, rmse, r, int(n)

# ========================== REALTIME TRACKER ==========================
class RealtimeTracker:
    def __init__(self, fs, cfg):
        self.fs, self.N, self.hop = float(fs), int(cfg["N"]), int(cfg["hop"])
        self.cfg = cfg
        
        # Janela Flattop (Coeficientes exatos mantidos)
        n = np.arange(self.N)
        self.win = (0.21557895 - 0.41663158*np.cos(2*np.pi*n/self.N) + 
                    0.277263158*np.cos(4*np.pi*n/self.N) - 0.083578947*np.cos(6*np.pi*n/self.N) + 
                    0.006947368*np.cos(8*np.pi*n/self.N))
        
        self.k_min = max(1, int(np.ceil(F_MIN_HZ * self.N / self.fs)) - 1)
        self.f_bins = np.fft.rfftfreq(self.N, d=1.0/self.fs)
        
        # Prior & EMA setup
        self.prior_min, self.prior_max = cfg["PRIOR_MIN_BPM"]/60.0, cfg["PRIOR_MAX_BPM"]/60.0
        self.f_prev = np.clip(cfg["PRIOR_INIT_BPM"]/60.0, self.prior_min, self.prior_max)
        self.ema_alpha = 1.0 - np.exp(-(self.hop/self.fs) / cfg["RT_EMA_TAU_SEC"]) if cfg["RT_EMA_TAU_SEC"] > 0 else 1.0
        self.ema_state = None
        
        self.buf, self.tbuf = deque(maxlen=self.N), deque(maxlen=self.N)
        self.t_out, self.hr_out = [], []
        self._cnt = 0

    def push(self, t, x):
        self.buf.append(x); self.tbuf.append(t)
        self._cnt += 1
        if len(self.buf) < self.N or self._cnt < self.hop: return
        self._cnt = 0

        # DSP
        xw = (np.array(self.buf) * self.win); xw -= np.mean(xw)
        P = np.abs(np.fft.rfft(xw, n=self.N))**2
        
        # Peak Tracking com Prior
        k_srch_min = max(self.k_min, int(np.ceil(self.prior_min * self.N / self.fs)))
        k_srch_max = min(len(P)-1, int(np.floor(self.prior_max * self.N / self.fs)))
        
        if k_srch_max > k_srch_min:
            W = np.exp(-0.5 * ((self.f_bins - self.f_prev) / self.cfg["PRIOR_SIGMA_HZ"])**2)
            k0 = k_srch_min + np.argmax((P * W)[k_srch_min : k_srch_max + 1])
        else:
            k0 = self.k_min + np.argmax(P[self.k_min:])

        f_est = (k0 + parabolic_peak(P, k0)) * self.fs / self.N
        
        # Post-process
        hr_inst = 60.0 * f_est
        self.ema_state = hr_inst if self.ema_state is None else (1-self.ema_alpha)*self.ema_state + self.ema_alpha*hr_inst
        
        # Update State & Lists
        self.f_prev = (1 - self.cfg["PRIOR_ALPHA"]) * self.f_prev + self.cfg["PRIOR_ALPHA"] * np.clip(f_est, self.prior_min, self.prior_max)
        self.t_out.append(self.tbuf[0] + 0.5*(self.N-1)/self.fs)
        self.hr_out.append(self.ema_state)

# ========================== MAIN FLOW ==========================
def main():
    sos, zi = design_sos_cheby2(FS_TARGET, CFG_REALTIME["FILT_ORD"], CFG_REALTIME["FILT_WI"], CFG_REALTIME["FILT_WF"], CFG_REALTIME["FILT_RS"])
    if sos is None: raise SystemExit("Erro Filtro")

    for ds in DATASETS:
        # 1. Leitura
        t_ms, xraw = read_radar_pandas(ds["radar"])
        t_pol, hr_pol = read_polar_pandas(ds["polar"])
        if len(t_ms) < 512 or len(t_pol) < 5:
            print(f"[REALTIME {ds['name']}] (dados insuficientes)")
            continue

        # 2. Pré-proc
        t_sec = (t_ms - t_ms[0]) / 1000.0
        t_u, x_u = segment_and_resample(t_sec, xraw - np.mean(xraw), FS_TARGET, GAP_THR)
        
        # 3. Simulação Realtime
        trk = RealtimeTracker(FS_TARGET, CFG_REALTIME)
        zi_rt = zi * x_u[0] if len(x_u) > 0 else zi
        for i in range(len(x_u)):
            y, zi_rt = sosfilt(sos, [x_u[i]], zi=zi_rt)
            trk.push(t_u[i], y[0])

        # 4. Métricas
        mae, rmse, r, n = calc_metrics(np.array(trk.hr_out), hr_pol, np.array(trk.t_out), t_pol)
        print(f"[REALTIME {ds['name']}] corr={r:.4f} | MAE={mae:.3f} | RMSE={rmse:.3f} | n={n}")

if __name__ == "__main__":
    main()