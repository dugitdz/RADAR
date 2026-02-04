import asyncio
import struct
import time
import threading
import csv
from collections import deque
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from bleak import BleakScanner, BleakClient

try:
    from scipy.signal import cheby2, filtfilt
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


DEVICE_NAME = "ESP32C3_IC"
CHAR_UUID   = "00001234-0000-1000-8000-00805f9b34fb"

CFG = dict(
    N=64,
    hop=32,
    f_min_hz=0.05,
    harm_on=1,
    harm_max=3,
    harm_tol=1,
    harm_ratio=0.20,
)

BP = dict(
    on=1,
    wi=0.80,
    wf=2.00,
    ford=4,
    Rs=30.0,
)

GUI_REFRESH_HZ  = 10.0
MAXP            = 1500
DEBUG_EVERY_N   = 50

CSV_ON          = 1
CSV_FLUSH_EVERY = 8
CSV_MAX_QUEUE   = 200000


def read_float_le(b4: bytes) -> float:
    u = (b4[0] | (b4[1] << 8) | (b4[2] << 16) | (b4[3] << 24)) & 0xFFFFFFFF
    return struct.unpack("<f", u.to_bytes(4, "little"))[0]


def parse_pkt_19(data: bytearray):
    if len(data) != 19:
        return None, f"len={len(data)}"
    if data[0] != 0x13:
        return None, f"hdr=0x{data[0]:02X}"
    total  = read_float_le(data[2:6])
    breath = read_float_le(data[6:10])
    heart  = read_float_le(data[10:14])
    t_ms   = read_float_le(data[15:19])
    if not (np.isfinite(total) and np.isfinite(breath) and np.isfinite(heart) and np.isfinite(t_ms)):
        return None, "nan/inf"
    return (t_ms, total, breath, heart), None


def flattopwin_periodic(N):
    a0, a1, a2, a3, a4 = 0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368
    n = np.arange(N, dtype=float)
    return (a0
            - a1*np.cos(2*np.pi*n/N)
            + a2*np.cos(4*np.pi*n/N)
            - a3*np.cos(6*np.pi*n/N)
            + a4*np.cos(8*np.pi*n/N))


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


def estimate_hr_window(x, fs, f_min_hz, harm_on, harm_max, harm_tol, harm_ratio):
    x = np.asarray(x, dtype=float).reshape(-1)
    N = x.size
    if N < 8 or (not np.isfinite(fs)) or fs <= 0:
        return np.nan
    w = flattopwin_periodic(N)
    xw = (x - np.mean(x)) * w
    X = np.fft.rfft(xw, n=N)
    P = np.abs(X)**2
    k_min_mat = int(np.ceil(f_min_hz * N / fs))
    k_min_mat = max(2, k_min_mat)
    k_min = max(1, k_min_mat - 1)
    if k_min >= P.size:
        return np.nan
    k0 = k_min + int(np.argmax(P[k_min:]))
    if harm_on:
        k0 = harmonic_guard_pick_subharmonic(P, k0, k_min, harm_max, harm_tol, harm_ratio)
    delta = parabolic_peak_interp_logP(P, k0)
    f_est = (k0 + delta) * fs / N
    return 60.0 * f_est


def _estimate_fs_from_t(tw):
    tw = np.asarray(tw, dtype=float).reshape(-1)
    if tw.size < 4:
        return np.nan
    dt = np.diff(tw)
    ok = np.isfinite(dt) & (dt > 0)
    if np.count_nonzero(ok) < max(3, tw.size // 4):
        return np.nan
    m = np.mean(dt[ok])
    if not np.isfinite(m) or m <= 0:
        return np.nan
    return 1.0 / m


def _try_filter_cheby2_filtfilt(x, fs, wi, wf, ford, Rs):
    if (not _HAVE_SCIPY) or (not np.isfinite(fs)) or fs <= 0:
        return None, "no_scipy_or_bad_fs"
    if (not np.isfinite(wi)) or (not np.isfinite(wf)) or wi <= 0 or wf <= 0 or wf <= wi:
        return None, "bad_wi_wf"
    nyq = 0.5 * fs
    if not np.isfinite(nyq) or nyq <= 0:
        return None, "bad_nyq"
    if wf >= nyq:
        return None, "wf>=nyq"
    wn = np.array([wi / nyq, wf / nyq], dtype=float)
    if np.any(~np.isfinite(wn)) or np.any(wn <= 0) or np.any(wn >= 1) or wn[1] <= wn[0]:
        return None, "bad_wn"
    if (ford % 2) != 0 or ford < 2:
        return None, "bad_ford"
    n = int(ford // 2)
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size < 8:
        return None, "too_short"
    try:
        b, a = cheby2(N=n, rs=float(Rs), Wn=wn, btype="bandpass", analog=False, output="ba")
    except Exception:
        return None, "cheby2_fail"
    try:
        padlen_default = 3 * (max(len(a), len(b)) - 1)
        if x.size <= padlen_default:
            return None, f"padlen_too_big({padlen_default})"
        y = filtfilt(b, a, x)
    except Exception:
        return None, "filtfilt_fail"
    if y.shape != x.shape or np.any(~np.isfinite(y)):
        return None, "nan_after_filter"
    return y, None


class State:
    def __init__(self):
        self.t0_ms = None
        self.last_t_ms = None

        self.t = deque(maxlen=200000)
        self.x = deque(maxlen=200000)

        self.hr_t = deque(maxlen=MAXP)
        self.hr   = deque(maxlen=MAXP)

        self.fs_t = deque(maxlen=MAXP)
        self.fs   = deque(maxlen=MAXP)

        self.sig_t    = deque(maxlen=8000)
        self.sig_raw  = deque(maxlen=8000)
        self.sig_filt = deque(maxlen=8000)

        self.good = 0
        self.bad  = 0
        self.ts_back = 0

        self.lock = threading.Lock()

        self.win_id = 0
        self.filter_ok = 0
        self.filter_fail = 0
        self.filter_last_err = ""

        self.csv_path = None
        self.csv_queue = deque(maxlen=CSV_MAX_QUEUE)
        self.csv_lock = threading.Lock()
        self.csv_flush_counter = 0


def _csv_init(st: State):
    if not CSV_ON:
        return
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    st.csv_path = f"ble_filter_debug_{ts}.csv"
    with open(st.csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["win_id", "t_center_s", "fs_win_hz", "bp_on", "bp_ok", "bp_err",
                    "t_s", "x_raw", "x_filt"])


def _csv_enqueue_window(st: State, win_id, t_center, fs, bp_on, bp_ok, bp_err, tw, x_raw, x_filt):
    if not CSV_ON or st.csv_path is None:
        return
    tw = np.asarray(tw, dtype=float).reshape(-1)
    x_raw = np.asarray(x_raw, dtype=float).reshape(-1)
    x_filt = np.asarray(x_filt, dtype=float).reshape(-1)
    n = min(tw.size, x_raw.size, x_filt.size)
    if n <= 0:
        return
    rows = []
    for i in range(n):
        rows.append([int(win_id), float(t_center), float(fs), int(bp_on), int(bp_ok), str(bp_err),
                     float(tw[i]), float(x_raw[i]), float(x_filt[i])])
    with st.csv_lock:
        for r in rows:
            st.csv_queue.append(r)
        st.csv_flush_counter += 1


def _csv_flush(st: State, force=False):
    if not CSV_ON or st.csv_path is None:
        return
    with st.csv_lock:
        if (not force) and (st.csv_flush_counter < CSV_FLUSH_EVERY):
            return
        if len(st.csv_queue) == 0:
            st.csv_flush_counter = 0
            return
        batch = list(st.csv_queue)
        st.csv_queue.clear()
        st.csv_flush_counter = 0
    try:
        with open(st.csv_path, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerows(batch)
    except Exception:
        pass


def process_if_possible(st: State, max_windows=2):
    N, hop = int(CFG["N"]), int(CFG["hop"])
    done = 0

    while done < int(max_windows):
        with st.lock:
            n = len(st.x)
            if n < N:
                return
            x_arr = np.asarray(st.x, dtype=float)
            t_arr = np.asarray(st.t, dtype=float)

        idx_end = n
        idx_start = idx_end - N
        if idx_start < 0:
            return

        tw = t_arr[idx_start:idx_end]
        xw_raw = x_arr[idx_start:idx_end]

        fs = _estimate_fs_from_t(tw)
        if not np.isfinite(fs) or fs <= 0:
            return

        t_center = float(tw[0] + 0.5 * (tw[-1] - tw[0]))

        bp_ok = 0
        bp_err = ""
        x_use = xw_raw

        if BP.get("on", 1):
            y, err = _try_filter_cheby2_filtfilt(
                xw_raw, fs, float(BP["wi"]), float(BP["wf"]), int(BP["ford"]), float(BP["Rs"])
            )
            if y is not None:
                x_use = y
                bp_ok = 1
            else:
                bp_err = err if err else "bp_unknown"

        t0 = time.perf_counter()
        hr = estimate_hr_window(
            x_use, fs, CFG["f_min_hz"],
            CFG["harm_on"], CFG["harm_max"], CFG["harm_tol"], CFG["harm_ratio"]
        )
        t1 = time.perf_counter()
        _proc_ms = 1000.0 * (t1 - t0)

        with st.lock:
            st.win_id += 1
            win_id = st.win_id
            if BP.get("on", 1):
                if bp_ok:
                    st.filter_ok += 1
                    st.filter_last_err = ""
                else:
                    st.filter_fail += 1
                    st.filter_last_err = bp_err

            st.fs_t.append(t_center); st.fs.append(float(fs))
            if np.isfinite(hr):
                st.hr_t.append(t_center); st.hr.append(float(hr))

            t_tail = tw[-hop:] if hop <= tw.size else tw
            raw_tail = xw_raw[-hop:] if hop <= xw_raw.size else xw_raw
            if bp_ok:
                filt_tail = x_use[-hop:] if hop <= x_use.size else x_use
            else:
                filt_tail = raw_tail

            for i in range(len(t_tail)):
                st.sig_t.append(float(t_tail[i]))
                st.sig_raw.append(float(raw_tail[i]))
                st.sig_filt.append(float(filt_tail[i]))

        x_filt_for_csv = x_use if bp_ok else xw_raw
        _csv_enqueue_window(
            st, win_id, t_center, fs, BP.get("on", 1), bp_ok, bp_err, tw, xw_raw, x_filt_for_csv
        )
        _csv_flush(st, force=False)

        if not np.isfinite(hr):
            return

        with st.lock:
            for _ in range(hop):
                if st.x:
                    st.x.popleft()
                    st.t.popleft()

        done += 1


def make_gui():
    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(11, 7), sharex=True)

    l_hr, = ax1.plot([], [], lw=1.5)
    p_hr, = ax1.plot([], [], linestyle="none", marker="o", markersize=4)
    ax1.set_ylabel("HR (bpm)")
    ax1.grid(True)

    l_fs, = ax2.plot([], [], lw=1.0)
    p_fs, = ax2.plot([], [], linestyle="none", marker="o", markersize=4)
    ax2.set_ylabel("Fs_win (Hz)")
    ax2.grid(True)

    l_raw,  = ax3.plot([], [], lw=1.0)
    l_filt, = ax3.plot([], [], lw=1.0)
    ax3.set_ylabel("heart_phase (raw/filt)")
    ax3.set_xlabel("t (s)")
    ax3.grid(True)
    ax3.legend(["raw", "filt"], loc="upper right")

    fig.tight_layout()
    return fig, (ax1, ax2, ax3), (l_hr, p_hr, l_fs, p_fs, l_raw, l_filt)


def gui_loop(st: State):
    fig, (ax1, ax2, ax3), (l_hr, p_hr, l_fs, p_fs, l_raw, l_filt) = make_gui()
    period = 1.0 / GUI_REFRESH_HZ

    while plt.fignum_exists(fig.number):
        with st.lock:
            t_hr = np.asarray(st.hr_t, dtype=float)
            y_hr = np.asarray(st.hr, dtype=float)

            t_fs = np.asarray(st.fs_t, dtype=float)
            y_fs = np.asarray(st.fs, dtype=float)

            ts = np.asarray(st.sig_t, dtype=float)
            xr = np.asarray(st.sig_raw, dtype=float)
            xf = np.asarray(st.sig_filt, dtype=float)

            good, bad, back = st.good, st.bad, st.ts_back
            fok, ffl, ferr = st.filter_ok, st.filter_fail, st.filter_last_err

        if t_hr.size > 0:
            l_hr.set_data(t_hr, y_hr)
            p_hr.set_data(t_hr, y_hr)
            ax1.relim(); ax1.autoscale_view()

        if t_fs.size > 0:
            l_fs.set_data(t_fs, y_fs)
            p_fs.set_data(t_fs, y_fs)
            ax2.relim(); ax2.autoscale_view()

        if ts.size > 1:
            l_raw.set_data(ts, xr)
            l_filt.set_data(ts, xf)
            ax3.relim(); ax3.autoscale_view()

        bp_txt = f"BP={'ON' if BP.get('on',1) else 'OFF'} {BP['wi']:.2f}-{BP['wf']:.2f}Hz Rs={BP['Rs']:.0f} ford={BP['ford']}"
        filt_txt = f"bp_ok={fok} bp_fail={ffl}" + (f" last={ferr}" if (ffl > 0 and ferr) else "")
        csv_txt = f"csv={st.csv_path}" if (CSV_ON and st.csv_path) else "csv=OFF"

        ax1.set_title(
            f"good={good} bad={bad} ts_back={back} | N={CFG['N']} hop={CFG['hop']} | {bp_txt} | {filt_txt} | {csv_txt}"
        )

        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(period)


def proc_loop(st: State):
    while True:
        try:
            process_if_possible(st, max_windows=2)
        except Exception:
            pass
        time.sleep(0.005)


async def main():
    st = State()
    _csv_init(st)

    gui_thread = threading.Thread(target=gui_loop, args=(st,), daemon=True)
    gui_thread.start()

    proc_thread = threading.Thread(target=proc_loop, args=(st,), daemon=True)
    proc_thread.start()

    print("Scanning BLE...")
    devices = await BleakScanner.discover(timeout=5.0)

    target = None
    for d in devices:
        if d.name == DEVICE_NAME:
            target = d
            break

    if not target:
        print(f"Não achei '{DEVICE_NAME}'. Vi:")
        for d in devices:
            print(" ", d.address, d.name)
        return

    def on_notify(sender: int, data: bytearray):
        rec, err = parse_pkt_19(data)
        if rec is None:
            st.bad += 1
            if st.bad % 20 == 0:
                try:
                    hx = data.hex()
                except Exception:
                    hx = "<nohex>"
                print(f"[BAD] {err} hex={hx}")
            return

        t_ms, total, breath, heart = rec

        if st.t0_ms is None:
            st.t0_ms = float(t_ms)

        if st.last_t_ms is not None and t_ms < st.last_t_ms:
            st.ts_back += 1
        st.last_t_ms = float(t_ms)

        t_sec = (float(t_ms) - float(st.t0_ms)) / 1000.0

        with st.lock:
            st.t.append(t_sec)
            st.x.append(float(heart))
            st.good += 1
            good = st.good

        if good % DEBUG_EVERY_N == 0:
            print(f"[OK] t_ms={t_ms:9.2f} t={t_sec:8.3f}s | total={total:+.3f} breath={breath:+.3f} heart={heart:+.3f}")

    print("Conectando em:", target.address, target.name)
    async with BleakClient(target.address) as client:
        if not client.is_connected:
            print("Falhou conectar.")
            return

        await client.start_notify(CHAR_UUID, on_notify)
        print("Recebendo (Ctrl+C pra sair)...")

        try:
            while True:
                await asyncio.sleep(1.0)
        except KeyboardInterrupt:
            pass

        await client.stop_notify(CHAR_UUID)
        _csv_flush(st, force=True)
        print("Parou.")


if __name__ == "__main__":
    asyncio.run(main())
