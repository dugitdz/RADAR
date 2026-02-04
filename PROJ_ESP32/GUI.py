import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime
from fractions import Fraction

from scipy.signal import resample_poly, get_window, butter, filtfilt, spectrogram
from scipy.interpolate import PchipInterpolator


# ===================== CONFIG (EDITÁVEL) ===================== #
PHASES_PATH   = r"C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases.csv"
POLAR_HR_PATH = r"C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH2.txt"

# Segmentação por buracos
GAP_THR_SEC = 0.5

# alvo de amostragem uniforme
FS_TARGET = 25.0  # Hz

# ===== rFFT por hop =====
NPERSEG     = 128
HOP_SAMPLES = NPERSEG // 5     # 80% overlap
WINDOW      = "hann"

# ===== ZERO-PHASE BANDPASS (filtfilt) =====
BP_ORDER = 2

# Heart band (Hz)
HEART_WI_HZ = 0.8
HEART_WF_HZ = 3.5

# Breath band (Hz)
BREATH_WI_HZ = 0.10
BREATH_WF_HZ = 0.60

# Remover DC / muito baixa freq (pra não travar em 0 Hz)
F_MIN_HZ = 0.05  # Hz

# ===== SUAVIZAÇÃO (em tempo no grid comum) =====
SMOOTH_HR_SEC = 5.0
SMOOTH_BR_SEC = 10.0

# ===== GRID COMUM (para interp/plot e movmean em tempo) =====
DT_GRID = 0.1  # s

# Comparação com Polar (métricas + overlay no espectro)
DO_POLAR_COMPARE = True

# eps pra log(P) na interp parabólica
EPS_PWR = 1e-20


# ===================== LEITURAS ===================== #
def read_phases_csv_matlab_style(path):
    """
    tpuro  = col 1 (ms)
    heart  = col 4 (index 3)
    breath = col 3 (index 2)
    """
    data = np.genfromtxt(path, delimiter=",", dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 4:
        raise ValueError(f"phases.csv precisa ter >=4 colunas, mas tem {data.shape[1]}")
    tpuro_ms = data[:, 0]
    heart    = data[:, 1]
    breath   = data[:, 2]
    return tpuro_ms, heart, breath


def read_txt_polar(path):
    df = pd.read_csv(path, sep=";", header=0, usecols=[0, 1])
    ts_str = df.iloc[:, 0].astype(str).to_list()
    hr = df.iloc[:, 1].to_numpy(dtype=np.float64)

    t_dt = []
    for s in ts_str:
        s = s.strip()
        try:
            t_dt.append(datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%f"))
        except ValueError:
            t_dt.append(datetime.strptime(s, "%Y-%m-%d %H:%M:%S.%f"))

    t0 = t_dt[0]
    t_sec = np.array([(x - t0).total_seconds() for x in t_dt], dtype=np.float64)

    # unique estável
    _, ia = np.unique(t_sec, return_index=True)
    ia = np.sort(ia)
    return t_sec[ia], hr[ia]


# ===================== UTIL ===================== #
def segment_by_gaps(t, gap_thr):
    dt = np.diff(t)
    breaks = np.where(dt > gap_thr)[0]
    n = len(t)
    if breaks.size == 0:
        return [(0, n - 1)]
    starts = np.concatenate(([0], breaks + 1))
    ends   = np.concatenate((breaks, [n - 1]))
    return list(zip(starts, ends))


def _estimate_fs_native(t_seg):
    dt = np.diff(t_seg)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    if dt.size < 3:
        return None
    fs = 1.0 / np.median(dt)
    if not np.isfinite(fs) or fs <= 0:
        return None
    return fs


def resample_poly_from_irregular(t_seg, x_seg, fs_target, fs_native_max_den=1000):
    """
    Irregular -> (uniform em fs_native via PCHIP) -> resample_poly -> fs_target.
    """
    t_seg = np.asarray(t_seg, dtype=np.float64)
    x_seg = np.asarray(x_seg, dtype=np.float64)

    if t_seg.size < 6:
        return None, None
    if np.any(np.diff(t_seg) <= 0):
        return None, None

    fs_native = _estimate_fs_native(t_seg)
    if fs_native is None:
        return None, None

    t0 = t_seg[0]
    t1 = t_seg[-1]
    if t1 <= t0:
        return None, None

    dt_native = 1.0 / fs_native
    tu_native = np.arange(t0, t1 + 1e-12, dt_native, dtype=np.float64)
    if tu_native.size < 8:
        return None, None

    try:
        f = PchipInterpolator(t_seg, x_seg, extrapolate=False)
        xu_native = f(tu_native)
    except Exception:
        return None, None

    ok = np.isfinite(xu_native)
    if np.sum(ok) < 8:
        return None, None
    tu_native = tu_native[ok]
    xu_native = xu_native[ok]

    xu_native = xu_native - np.mean(xu_native)

    ratio = Fraction(fs_target / fs_native).limit_denominator(fs_native_max_den)
    up, down = ratio.numerator, ratio.denominator
    if up <= 0 or down <= 0:
        return None, None

    xu_target = resample_poly(xu_native, up, down).astype(np.float64)
    tu_target = t0 + np.arange(xu_target.size, dtype=np.float64) / float(fs_target)

    if xu_target.size < 8:
        return None, None

    return tu_target, xu_target


def zero_phase_bandpass(x, fs, wi_hz, wf_hz, order):
    x = np.asarray(x, dtype=np.float64)
    if x.size < 8:
        return x

    nyq = 0.5 * fs
    wi = float(wi_hz) / nyq
    wf = float(wf_hz) / nyq

    wi = max(1e-6, min(wi, 0.999999))
    wf = max(1e-6, min(wf, 0.999999))
    if wf <= wi:
        return x

    b, a = butter(int(order), [wi, wf], btype="bandpass")

    padlen = 3 * (max(len(a), len(b)) - 1)
    if x.size <= padlen:
        return x

    try:
        return filtfilt(b, a, x)
    except Exception:
        return x


def parabolic_peak_interp_logP(P, k):
    if k <= 0 or k >= (P.size - 1):
        return 0.0

    a0 = float(P[k - 1])
    b0 = float(P[k])
    c0 = float(P[k + 1])

    la = np.log(a0 + EPS_PWR)
    lb = np.log(b0 + EPS_PWR)
    lc = np.log(c0 + EPS_PWR)

    denom = (la - 2.0 * lb + lc)
    if denom == 0 or not np.isfinite(denom):
        return 0.0

    delta = 0.5 * (la - lc) / denom
    if not np.isfinite(delta):
        delta = 0.0
    return float(np.clip(delta, -0.75, 0.75))


def rfft_fundamental_track_hop(tu, xu, fs, nperseg, hop_samples, window, f_min_hz):
    xu = np.asarray(xu, dtype=np.float64)
    tu = np.asarray(tu, dtype=np.float64)

    if xu.size < nperseg or tu.size != xu.size:
        return None, None

    hop_samples = int(max(1, hop_samples))
    nperseg = int(nperseg)

    win = get_window(window, nperseg, fftbins=True).astype(np.float64)
    nfft = nperseg

    k_min = int(np.ceil(float(f_min_hz) * nfft / float(fs)))
    k_min = max(1, k_min)

    t_frames = []
    fpk = []

    for i0 in range(0, xu.size - nperseg + 1, hop_samples):
        xw = xu[i0:i0+nperseg] * win
        X = np.fft.rfft(xw, n=nfft)
        P = (np.abs(X) ** 2).astype(np.float64)

        if k_min >= P.size:
            continue

        k0 = k_min + int(np.argmax(P[k_min:]))
        delta = parabolic_peak_interp_logP(P, k0)
        f_est = (k0 + delta) * fs / nfft

        t_center = tu[i0] + 0.5 * (nperseg - 1) / fs
        t_frames.append(t_center)
        fpk.append(f_est)

    if len(t_frames) == 0:
        return None, None

    return np.array(t_frames, dtype=np.float64), np.array(fpk, dtype=np.float64)


def interp_common(t_src, y_src, t_dst):
    t_src = np.asarray(t_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    ok = np.isfinite(t_src) & np.isfinite(y_src)
    out = np.full_like(t_dst, np.nan, dtype=np.float64)
    if np.sum(ok) < 2:
        return out
    out[:] = np.interp(t_dst, t_src[ok], y_src[ok], left=np.nan, right=np.nan)
    return out


def movmean_nan_uniform(y, win_n):
    y = np.asarray(y, dtype=np.float64)
    if win_n <= 1:
        return y.copy()

    valid = np.isfinite(y).astype(np.float64)
    y0 = np.where(np.isfinite(y), y, 0.0)

    k = np.ones(int(win_n), dtype=np.float64)
    num = np.convolve(y0, k, mode="same")
    den = np.convolve(valid, k, mode="same")

    out = np.full_like(y, np.nan, dtype=np.float64)
    m = den > 0
    out[m] = num[m] / den[m]
    return out


def smooth_in_time_to_grid(t_src, y_src, t_grid, smooth_sec):
    y_i = interp_common(t_src, y_src, t_grid)
    win_n = int(max(1, round(float(smooth_sec) / float(DT_GRID))))
    return movmean_nan_uniform(y_i, win_n)


def metrics(yhat, yref):
    m = np.isfinite(yhat) & np.isfinite(yref)
    if np.sum(m) < 3:
        return np.nan, np.nan, np.nan, 0
    e = yhat[m] - yref[m]
    mae  = float(np.mean(np.abs(e)))
    rmse = float(np.sqrt(np.mean(e**2)))
    corr = float(np.corrcoef(yhat[m], yref[m])[0, 1])
    return mae, rmse, corr, int(np.sum(m))


# ===================== MAIN ===================== #
def main():
    tpuro_ms, heart, breath = read_phases_csv_matlab_style(PHASES_PATH)
    t = (tpuro_ms - tpuro_ms[0]) / 1000.0
    segments = segment_by_gaps(t, GAP_THR_SEC)

    # frames (HR/BR)
    tH_all, fH_all = [], []
    tB_all, fB_all = [], []

    # maior segmento HEART filtrado (pra espectrograma)
    best_spec_len = -1
    tu_spec = None
    hu_spec_f = None

    for s0, e0 in segments:
        tseg = t[s0:e0+1]
        hseg = heart[s0:e0+1]
        bseg = breath[s0:e0+1]

        # HEART
        tu, hu = resample_poly_from_irregular(tseg, hseg, FS_TARGET)
        if tu is not None and tu.size >= NPERSEG:
            hu_f = zero_phase_bandpass(hu, FS_TARGET, HEART_WI_HZ, HEART_WF_HZ, BP_ORDER)

            if hu_f.size > best_spec_len:
                best_spec_len = hu_f.size
                tu_spec = tu.copy()
                hu_spec_f = hu_f.copy()

            t_frames, fpk = rfft_fundamental_track_hop(
                tu, hu_f, FS_TARGET, NPERSEG, HOP_SAMPLES, WINDOW, F_MIN_HZ
            )
            if t_frames is not None:
                tH_all.append(t_frames)
                fH_all.append(fpk)

        # BREATH
        tu, bu = resample_poly_from_irregular(tseg, bseg, FS_TARGET)
        if tu is not None and tu.size >= NPERSEG:
            bu_f = zero_phase_bandpass(bu, FS_TARGET, BREATH_WI_HZ, BREATH_WF_HZ, BP_ORDER)
            t_frames, fpk = rfft_fundamental_track_hop(
                tu, bu_f, FS_TARGET, NPERSEG, HOP_SAMPLES, WINDOW, F_MIN_HZ
            )
            if t_frames is not None:
                tB_all.append(t_frames)
                fB_all.append(fpk)

    if len(tH_all) == 0:
        raise RuntimeError("Não consegui extrair HR: segmentos muito curtos ou t inválido.")

    # concat HR
    t_hr_frames = np.concatenate(tH_all)
    hr_bpm_frames = 60.0 * np.concatenate(fH_all)

    # concat BR
    if len(tB_all) > 0:
        t_br_frames = np.concatenate(tB_all)
        br_brpm_frames = 60.0 * np.concatenate(fB_all)
    else:
        t_br_frames = np.array([], dtype=np.float64)
        br_brpm_frames = np.array([], dtype=np.float64)

    # ===================== POLAR ===================== #
    t_polar, hr_polar = None, None
    if DO_POLAR_COMPARE:
        try:
            t_polar, hr_polar = read_txt_polar(POLAR_HR_PATH)
        except Exception:
            t_polar, hr_polar = None, None

    # ===================== GRID comum (HR suavizado) ===================== #
    t_common = np.arange(np.nanmin(t_hr_frames), np.nanmax(t_hr_frames) + 1e-12, DT_GRID, dtype=np.float64)
    hr_radar_s = smooth_in_time_to_grid(t_hr_frames, hr_bpm_frames, t_common, SMOOTH_HR_SEC)

    if t_polar is not None and hr_polar is not None:
        hr_polar_i = interp_common(t_polar, hr_polar, t_common)
        mae, rmse, corr, n = metrics(hr_radar_s, hr_polar_i)
        print(f"Radar x Polar | N={n} | MAE={mae:.3f} bpm | RMSE={rmse:.3f} bpm | Corr={corr:.3f}")

    # ===================== PLOT HR (opcional, mas útil) ===================== #
    plt.figure()
    plt.plot(t_common, hr_radar_s, label="Radar HR (movmean)")
    if t_polar is not None and hr_polar is not None:
        plt.plot(t_polar, hr_polar, label="Polar HR")
    plt.xlabel("Tempo (s)")
    plt.ylabel("HR (bpm)")
    plt.title("HR (Radar vs Polar)")
    plt.grid(True)
    plt.legend()

    # ===================== PLOT BR (opcional) ===================== #
    if t_br_frames.size > 0:
        tbr_grid = np.arange(np.nanmin(t_br_frames), np.nanmax(t_br_frames) + 1e-12, DT_GRID, dtype=np.float64)
        br_radar_s = smooth_in_time_to_grid(t_br_frames, br_brpm_frames, tbr_grid, SMOOTH_BR_SEC)

        plt.figure()
        plt.plot(tbr_grid, br_radar_s, label="Radar BR (movmean)")
        plt.xlabel("Tempo (s)")
        plt.ylabel("RR (brpm)")
        plt.title("BR (Radar)")
        plt.grid(True)
        plt.legend()

    # ===================== ESPECTROGRAMA + POLAR POR CIMA (NO MESMO AX) ===================== #
    if tu_spec is not None and hu_spec_f is not None and hu_spec_f.size >= NPERSEG:
        noverlap = max(0, NPERSEG - HOP_SAMPLES)

        f, tt, Sxx = spectrogram(
            hu_spec_f,
            fs=FS_TARGET,
            window=get_window(WINDOW, NPERSEG),
            nperseg=NPERSEG,
            noverlap=noverlap,
            detrend=False,
            scaling="density",
            mode="psd",
        )

        tt_abs = tu_spec[0] + tt

        # 0..-40 dB relativo ao pico
        Sxx_db = 10.0 * np.log10(Sxx + 1e-30)
        Sxx_db = Sxx_db - np.nanmax(Sxx_db)
        Sxx_db = np.clip(Sxx_db, -40.0, 0.0)

        fig, ax = plt.subplots(figsize=(10, 5))
        pcm = ax.pcolormesh(tt_abs, f, Sxx_db, shading="auto", vmin=-40.0, vmax=0.0)

        ax.set_ylim(0.5, 3.0)  # ajuste se quiser
        ax.set_xlabel("Tempo (s)")
        ax.set_ylabel("Frequência (Hz)")
        ax.set_title(
            f"Espectrograma (HEART filtrado) | Fs={FS_TARGET:.1f}Hz | N={NPERSEG} | overlap={noverlap}\n"
            f"BP: {HEART_WI_HZ:.2f}–{HEART_WF_HZ:.2f} Hz | escala: 0..-40 dB (rel. ao pico)"
        )

        cb = fig.colorbar(pcm, ax=ax, pad=0.02)
        cb.set_label("PSD (dB rel. ao pico)")

        # Polar por cima no MESMO PLOT usando eixo Y secundário em bpm (alinhado Hz*60)
        ax_bpm = ax.twinx()
        y0_hz, y1_hz = ax.get_ylim()
        ax_bpm.set_ylim(y0_hz * 60.0, y1_hz * 60.0)
        ax_bpm.set_ylabel("HR (bpm)")

        if t_polar is not None and hr_polar is not None:
            mpol = (
                np.isfinite(t_polar) & np.isfinite(hr_polar) &
                (t_polar >= tt_abs.min()) & (t_polar <= tt_abs.max())
            )
            if np.any(mpol):
                ax_bpm.plot(t_polar[mpol], hr_polar[mpol], linewidth=2.0, label="Polar HR")
                ax_bpm.legend(loc="upper right")

        plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
