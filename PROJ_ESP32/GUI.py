import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime
from scipy.signal import cheby2, filtfilt
from scipy.interpolate import PchipInterpolator

# ========================== PATHS (3 DUPLAS) ==========================
BASE = r"C:\Users\eduar\UTFPR\IC\RADAR\PROJ_ESP32\output\\"
DATASETS = [
    {"name": "H2",  "radar": BASE + "phases.csv",     "polar": BASE + "POLARH2.txt"},
    {"name": "H3",  "radar": BASE + "phases_raw.csv", "polar": BASE + "POLARH3.txt"},
    {"name": "H10", "radar": BASE + "phase.csv",      "polar": BASE + "POLARH10.txt"},
]

# ===================== FIXOS =====================
GAP_THR_SEC = 0.6
FS_TARGET   = 25.0

COL_T_MS  = 0
COL_HEART = 1

DT_GRID       = 0.1
SMOOTH_HR_SEC = 7.0
F_MIN_HZ      = 0.05

MIN_POINTS_AFTER_RESAMPLE = 32
MIN_METRIC_SAMPLES        = 10

# ===================== CFG DEFINITIVA =====================
CFG = {
    "N": 32,
    "hop": 16,
    "WRAP": 0,

    "win": "flattop",   # fixo
    "ftype": "cheby2",  # fixo
    "ford": 4,          # par
    "wi": 0.80,
    "wf": 2.00,

    "harm_on": 1,
    "harm_max": 3,
    "harm_tol": 1,
    "harm_ratio": 0.20,
}

PLOT_ON = True
PLOT_NONBLOCK = True


# ========================== FUNÇÕES ==========================
def _extract_floats_from_line(line: str):
    # extrai floats estilo readmatrix (pega inclusive notação científica)
    # sem usar regex (pra ficar leve): troca separadores por espaço e split
    for ch in [",", ";", "\t", "|"]:
        line = line.replace(ch, " ")
    parts = line.strip().split()
    vals = []
    for p in parts:
        try:
            vals.append(float(p))
        except Exception:
            pass
    return vals

def read_radar_like_readmatrix(path):
    # emula readmatrix: ignora linhas ruins e pega as 2 primeiras colunas numéricas
    t_ms = []
    x = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            vals = _extract_floats_from_line(line)
            if len(vals) >= 2:
                t_ms.append(vals[0])
                x.append(vals[1])
    t_ms = np.asarray(t_ms, dtype=float)
    x = np.asarray(x, dtype=float)
    return t_ms, x

def interp1_index_linear_extrap(x):
    # MATLAB:
    # x(~ok)=interp1(find(ok),x(ok),find(~ok),'linear','extrap')
    x = np.asarray(x, dtype=float).reshape(-1)
    ok = np.isfinite(x)
    if np.all(ok):
        return x

    idx_ok = np.flatnonzero(ok)
    idx_q  = np.flatnonzero(~ok)
    y_ok   = x[idx_ok]

    if idx_ok.size == 0:
        return np.zeros_like(x)

    if idx_ok.size == 1:
        x[idx_q] = y_ok[0]
        return x

    # interp dentro do range
    x[idx_q] = np.interp(idx_q, idx_ok, y_ok)

    # extrap linear à esquerda
    left = idx_q[idx_q < idx_ok[0]]
    if left.size > 0:
        i0, i1 = idx_ok[0], idx_ok[1]
        y0, y1 = y_ok[0], y_ok[1]
        m = (y1 - y0) / (i1 - i0)
        x[left] = y0 + m * (left - i0)

    # extrap linear à direita
    right = idx_q[idx_q > idx_ok[-1]]
    if right.size > 0:
        i0, i1 = idx_ok[-2], idx_ok[-1]
        y0, y1 = y_ok[-2], y_ok[-1]
        m = (y1 - y0) / (i1 - i0)
        x[right] = y1 + m * (right - i1)

    return x

def prep_nowrap(x):
    x = np.asarray(x, dtype=float).reshape(-1)
    x = interp1_index_linear_extrap(x)
    return x - np.nanmean(x)

def wrap_phase(ph):
    ph = np.asarray(ph, dtype=float).reshape(-1)
    ph = interp1_index_linear_extrap(ph)
    phw = (ph + np.pi) % (2*np.pi) - np.pi
    return phw - np.nanmean(phw)

def segment_by_gaps(t, gap_thr):
    t = np.asarray(t, dtype=float).reshape(-1)
    n = t.size
    if n < 2:
        return np.array([[0, max(0, n-1)]], dtype=int)
    br = np.flatnonzero(np.diff(t) > gap_thr)
    if br.size == 0:
        return np.array([[0, n-1]], dtype=int)
    starts = np.r_[0, br + 1]
    ends   = np.r_[br, n-1]
    return np.c_[starts, ends].astype(int)

def resample_segments_irregular(t, x, segments, fs_target, method="linear"):
    t = np.asarray(t, dtype=float).reshape(-1)
    x = np.asarray(x, dtype=float).reshape(-1)
    dt = 1.0 / fs_target

    t_u_list, x_u_list = [], []

    for i0, i1 in segments:
        ts = t[i0:i1+1]
        xs = x[i0:i1+1]

        # unique(ts,'stable')
        _, iu = np.unique(ts, return_index=True)
        iu = np.sort(iu)
        ts = ts[iu]
        xs = xs[iu]

        if ts.size < 4:
            continue

        # MATLAB: (ts(1):dt:ts(end))'
        tnew = np.arange(ts[0], ts[-1] + 1e-12, dt)

        if method == "pchip":
            itp = PchipInterpolator(ts, xs, extrapolate=True)
            xnew = itp(tnew)
        else:
            xnew = np.interp(tnew, ts, xs)

        ok = np.isfinite(tnew) & np.isfinite(xnew)
        if np.any(ok):
            t_u_list.append(tnew[ok])
            x_u_list.append(xnew[ok])

    if not t_u_list:
        return np.array([]), np.array([])

    t_u = np.concatenate(t_u_list)
    x_u = np.concatenate(x_u_list)

    # unique(t_u,'stable')
    _, iu = np.unique(t_u, return_index=True)
    iu = np.sort(iu)
    return t_u[iu], x_u[iu]

def zero_phase_bandpass_cheby2(x, fs, wi_hz, wf_hz, order_final):
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size < 16:
        return x
    if wf_hz >= fs/2 or wi_hz <= 0 or wi_hz >= wf_hz or (order_final % 2) != 0:
        return x

    n = order_final // 2
    Rs = 30.0
    wn = np.array([wi_hz, wf_hz]) / (fs/2)

    b, a = cheby2(n, Rs, wn, btype="bandpass", output="ba")

    padlen = 3 * (max(len(a), len(b)) - 1)
    if x.size <= padlen:
        return x

    return filtfilt(b, a, x)

def flattopwin_periodic(N):
    # MATLAB flattopwin(N,'periodic') (coef padrão)
    # w[n] = a0 - a1*cos(2πn/N) + a2*cos(4πn/N) - a3*cos(6πn/N) + a4*cos(8πn/N)
    a0 = 0.21557895
    a1 = 0.41663158
    a2 = 0.277263158
    a3 = 0.083578947
    a4 = 0.006947368
    n = np.arange(N, dtype=float)
    w = (a0
         - a1*np.cos(2*np.pi*n/N)
         + a2*np.cos(4*np.pi*n/N)
         - a3*np.cos(6*np.pi*n/N)
         + a4*np.cos(8*np.pi*n/N))
    return w

def parabolic_peak_interp_logP(P, k):
    if k <= 0 or k >= (len(P) - 1):
        return 0.0
    la = np.log(P[k-1] + 1e-20)
    lb = np.log(P[k]   + 1e-20)
    lc = np.log(P[k+1] + 1e-20)
    den = (la - 2*lb + lc)
    if den == 0 or not np.isfinite(den):
        return 0.0
    d = 0.5 * (la - lc) / den
    if not np.isfinite(d):
        d = 0.0
    return float(np.clip(d, -0.75, 0.75))

def harmonic_guard_pick_subharmonic(P, k0, k_min, max_order, tol_bins, ratio_min):
    # FIEL AO MATLAB (k_est sempre baseado no k0 original)
    k_pick = int(k0)
    P0 = P[k_pick]
    if (not np.isfinite(P0)) or (P0 <= 0):
        return k_pick

    k0_base = int(k0)

    for m in range(2, int(max_order) + 1):
        k_est = int(np.round(k0_base / m))
        if k_est < k_min:
            continue

        lo = max(int(k_min), k_est - int(tol_bins))
        hi = min(len(P) - 1, k_est + int(tol_bins))
        if hi < lo:
            continue

        seg = P[lo:hi+1]
        idx = int(np.argmax(seg))
        ksub = lo + idx
        Psub = seg[idx]

        if np.isfinite(Psub) and (Psub >= ratio_min * P0):
            k_pick = ksub
            P0 = Psub

    return int(k_pick)

def rfft_fundamental_track_hop(tu, xu, fs, nperseg, hop, f_min_hz,
                               harm_on, harm_max, harm_tol, harm_ratio):
    tu = np.asarray(tu, dtype=float).reshape(-1)
    xu = np.asarray(xu, dtype=float).reshape(-1)

    if xu.size < nperseg or tu.size != xu.size:
        return np.array([]), np.array([])

    # MATLAB: flattopwin(nperseg,'periodic')
    win = flattopwin_periodic(nperseg)

    nfft = nperseg
    half = nfft//2 + 1

    # MATLAB:
    # k_min_mat = ceil(f_min*nfft/fs); k_min_mat=max(2,k_min_mat)
    # python index = k_mat-1
    k_min_mat = int(np.ceil(f_min_hz * nfft / fs))
    k_min_mat = max(2, k_min_mat)
    k_min = max(1, k_min_mat - 1)

    t_frames = []
    fpk_hz = []

    idx0 = 0
    while (idx0 + nperseg) <= xu.size:
        xw = xu[idx0:idx0+nperseg] * win
        xw = xw - np.nanmean(xw)

        X = np.fft.rfft(xw, n=nfft)
        P = np.abs(X[:half])**2

        if k_min >= P.size:
            idx0 += hop
            continue

        rel = int(np.argmax(P[k_min:]))
        k0 = k_min + rel

        if harm_on:
            k0 = harmonic_guard_pick_subharmonic(P, k0, k_min, harm_max, harm_tol, harm_ratio)

        delta = parabolic_peak_interp_logP(P, k0)

        # MATLAB: f = (k0_mat - 1 + delta)*fs/nfft
        # python: k0_py = k0_mat-1 => f = (k0_py + delta)*fs/nfft
        f_est = (k0 + delta) * fs / nfft

        t_center = tu[idx0] + 0.5*(nperseg-1)/fs

        t_frames.append(t_center)
        fpk_hz.append(f_est)

        idx0 += hop

    return np.asarray(t_frames), np.asarray(fpk_hz)

def interp_common(t_src, y_src, t_dst, kind="linear"):
    t_src = np.asarray(t_src, dtype=float).reshape(-1)
    y_src = np.asarray(y_src, dtype=float).reshape(-1)
    t_dst = np.asarray(t_dst, dtype=float).reshape(-1)

    ok = np.isfinite(t_src) & np.isfinite(y_src)
    out = np.full_like(t_dst, np.nan, dtype=float)
    if np.count_nonzero(ok) < 2:
        return out

    ts = t_src[ok]
    ys = y_src[ok]

    iu = np.argsort(ts)
    ts = ts[iu]
    ys = ys[iu]

    ts_u, ia = np.unique(ts, return_index=True)
    ys_u = ys[ia]
    if ts_u.size < 2:
        return out

    if kind == "pchip":
        itp = PchipInterpolator(ts_u, ys_u, extrapolate=False)
        out = itp(t_dst)
    else:
        out = np.interp(t_dst, ts_u, ys_u, left=np.nan, right=np.nan)

    return out

def smooth_in_time_to_grid(t_src, y_src, t_grid, dt_grid, smooth_sec):
    y_i = interp_common(t_src, y_src, t_grid, kind="linear")
    win_n = max(1, int(np.round(smooth_sec / dt_grid)))

    m = np.isfinite(y_i).astype(float)
    y0 = np.nan_to_num(y_i, nan=0.0)
    ker = np.ones(win_n, dtype=float)

    num = np.convolve(y0, ker, mode="same")
    den = np.convolve(m,  ker, mode="same")

    out = np.full_like(y_i, np.nan, dtype=float)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out

def metrics(yhat, yref):
    yhat = np.asarray(yhat, dtype=float).reshape(-1)
    yref = np.asarray(yref, dtype=float).reshape(-1)
    m = np.isfinite(yhat) & np.isfinite(yref)
    n = int(np.count_nonzero(m))
    if n < 3:
        return np.nan, np.nan, np.nan, n
    e = yhat[m] - yref[m]
    mae = float(np.mean(np.abs(e)))
    rmse = float(np.sqrt(np.mean(e*e)))

    yh = yhat[m]
    yr = yref[m]
    if np.std(yh) == 0 or np.std(yr) == 0:
        r = np.nan
    else:
        r = float(np.corrcoef(yh, yr)[0, 1])
    return mae, rmse, r, n

def _parse_dt(s):
    s = str(s).strip()
    fmts = [
        "%Y-%m-%dT%H:%M:%S.%f",
        "%Y-%m-%d %H:%M:%S.%f",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
    ]
    for f in fmts:
        try:
            return datetime.strptime(s, f)
        except Exception:
            pass
    try:
        return datetime.fromisoformat(s.replace("Z", ""))
    except Exception:
        raise ValueError(f"Timestamp Polar inválido: {s}")

def read_txt_polar(path):
    df = pd.read_csv(path, sep=";", header=0, usecols=[0, 1], engine="python")
    ts = df.iloc[:, 0].astype(str).to_numpy()
    hr = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)

    ok = np.isfinite(hr)
    ts = ts[ok]
    hr = hr[ok]
    if ts.size == 0:
        return np.array([]), np.array([])

    tdt = np.array([_parse_dt(s) for s in ts], dtype=object)
    t0 = tdt[0]
    t_sec = np.array([(tt - t0).total_seconds() for tt in tdt], dtype=float)

    _, ia = np.unique(t_sec, return_index=True)
    ia = np.sort(ia)
    return t_sec[ia], hr[ia]


# ===================== RUN =====================
print("\n==================== CONFIG DEFINITIVA ====================")
print(f"N={CFG['N']} hop={CFG['hop']} WRAP={CFG['WRAP']} win={CFG['win']} | "
      f"filt={CFG['ftype']} ord={CFG['ford']} wi={CFG['wi']:.2f} wf={CFG['wf']:.2f} | "
      f"harm={CFG['harm_on']} (max={CFG['harm_max']} tol={CFG['harm_tol']} ratio={CFG['harm_ratio']:.2f})")

Corrs, MAEs, RMSEs = [], [], []
used = 0

print("\n==================== RUN (DEFINITIVO) ====================")
for ds in DATASETS:
    print(f"\n[DATASET] {ds['name']}")

    # >>> aqui resolve o H3 (readmatrix-like)
    t_ms, xraw = read_radar_like_readmatrix(ds["radar"])
    if t_ms.size < 4:
        print(" [SKIP] Radar vazio/poucos pontos.")
        continue

    t0 = (t_ms - t_ms[0]) / 1000.0
    seg0 = segment_by_gaps(t0, GAP_THR_SEC)

    try:
        t_polar, HR_polar = read_txt_polar(ds["polar"])
        hasPolar = (t_polar.size > 0) and (HR_polar.size > 0)
    except Exception:
        hasPolar = False

    if not hasPolar:
        print(" [SKIP] Polar indisponível.")
        continue

    x0 = wrap_phase(xraw) if CFG["WRAP"] else prep_nowrap(xraw)

    t_u, x_u = resample_segments_irregular(t0, x0, seg0, FS_TARGET, method="linear")
    if x_u.size < max(MIN_POINTS_AFTER_RESAMPLE, CFG["N"]):
        print(f" [SKIP] Poucos pontos após resample: {x_u.size}")
        continue

    x_f = zero_phase_bandpass_cheby2(x_u, FS_TARGET, CFG["wi"], CFG["wf"], CFG["ford"])

    t_frames, f_hz = rfft_fundamental_track_hop(
        t_u, x_f, FS_TARGET,
        CFG["N"], CFG["hop"], F_MIN_HZ,
        CFG["harm_on"], CFG["harm_max"], CFG["harm_tol"], CFG["harm_ratio"]
    )

    if t_frames.size == 0:
        print(" [SKIP] Track vazio.")
        continue

    hr_bpm = 60.0 * f_hz

    t_start = max(np.min(t_frames), np.min(t_polar))
    t_end   = min(np.max(t_frames), np.max(t_polar))
    if not (np.isfinite(t_start) and np.isfinite(t_end) and (t_end > t_start)):
        print(" [SKIP] Janela comum inválida.")
        continue

    t_common = np.arange(t_start, t_end + 1e-12, DT_GRID)

    HR_radar_grid = smooth_in_time_to_grid(t_frames, hr_bpm, t_common, DT_GRID, SMOOTH_HR_SEC)
    HR_polar_grid = interp_common(t_polar, HR_polar, t_common, kind="linear")

    MAE, RMSE, R, n = metrics(HR_radar_grid, HR_polar_grid)
    if not (np.isfinite(MAE) and np.isfinite(R) and (n >= MIN_METRIC_SAMPLES)):
        print(f" [SKIP] Métrica inválida (n={n}).")
        continue

    used += 1
    Corrs.append(R)
    MAEs.append(MAE)
    RMSEs.append(RMSE)

    print(f" OK | corr={R:.4f} | MAE={MAE:.3f} | RMSE={RMSE:.3f} | n={n}")

    if PLOT_ON:
        plt.figure(figsize=(10, 4))
        plt.plot(t_common, HR_polar_grid, label="Polar (grid)")
        plt.plot(t_common, HR_radar_grid, label="Radar (grid)")
        plt.grid(True)
        plt.xlabel("t (s)")
        plt.ylabel("HR (bpm)")
        plt.title(f"DEFINITIVO {ds['name']} | corr={R:.3f} MAE={MAE:.2f} RMSE={RMSE:.2f}")
        plt.legend()
        plt.tight_layout()
        if PLOT_NONBLOCK:
            plt.show(block=False)
            plt.pause(0.01)
        else:
            plt.show()

print("\n==================== SUMMARY ====================")
if used == 0:
    print("Nenhum dataset válido.")
else:
    print(f"Datasets usados: {used}/{len(DATASETS)}")
    print(f"corr_mean={np.mean(Corrs):.4f} | MAE_mean={np.mean(MAEs):.3f} | RMSE_mean={np.mean(RMSEs):.3f}")
