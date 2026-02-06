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


DEVICE_NAME = "ESP32C3_IC"
CHAR_UUID   = "00001234-0000-1000-8000-00805f9b34fb"

# ===================== CONFIG =====================
CFG = dict(
    
    N=64,
    hop=32,

    # --- HR (bpm) ---
    HR_COL=3,          # 1=total, 2=breath, 3=heart
    hr_f_min_hz=0.05,
    hr_wi=0.60,
    hr_wf=2.50,
    hr_harm_on=1,
    hr_harm_max=3,
    hr_harm_tol=1,
    hr_harm_ratio=0.20,

    # --- BR (breaths per minute) ---
    BR_COL=2,          # 1=total, 2=breath, 3=heart
    br_f_min_hz=0.05,  
    br_wi=0.10,       
    br_wf=0.50,      
    br_harm_on=0,      
    br_harm_max=3,
    br_harm_tol=1,
    br_harm_ratio=0.25,
)

GUI_REFRESH_HZ  = 10.0
MAXP            = 1500
DEBUG_EVERY_N   = 50

CSV_ON          = 1
CSV_FLUSH_EVERY = 8
CSV_MAX_QUEUE   = 200000


# ===================== BLE PACKET =====================
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


# ===================== DSP =====================
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


def detrend_linear(x):
    x = np.asarray(x, dtype=float).reshape(-1)
    N = x.size
    if N < 3:
        return x
    n = np.arange(N, dtype=float)
    A = np.vstack([n, np.ones(N)]).T
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    trend = A @ coef
    return x - trend


def estimate_rate_window(x, fs, f_min_hz, harm_on, harm_max, harm_tol, harm_ratio, wi, wf):
    """
    Estima frequência fundamental no band [wi,wf] e retorna em "per-minute" (60*f).
    Serve tanto pra HR (bpm) quanto BR (brpm).
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    N = x.size
    if N < 8 or (not np.isfinite(fs)) or fs <= 0:
        return np.nan

    x = detrend_linear(x)

    w = flattopwin_periodic(N)
    xw = (x - np.mean(x)) * w

    X = np.fft.rfft(xw, n=N)
    P = np.abs(X)**2

    k_min_mat = int(np.ceil(f_min_hz * N / fs))
    k_min_mat = max(2, k_min_mat)
    k_min = max(1, k_min_mat - 1)

    if (not np.isfinite(wi)) or (not np.isfinite(wf)) or wi <= 0 or wf <= wi:
        return np.nan

    # CORRIGIDO: k_lo com CEIL (pra não colar em bin baixo demais)
    k_lo = int(np.ceil(wi * N / fs))
    k_hi = int(np.ceil(wf * N / fs))
    k_lo = max(k_min, k_lo)
    k_hi = min(P.size - 1, k_hi)

    if k_hi <= k_lo or k_lo >= P.size:
        return np.nan

    k0 = k_lo + int(np.argmax(P[k_lo:k_hi+1]))

    # reject mínimo (só se colar no limite inferior)
    if k0 <= k_lo:
        return np.nan

    if harm_on:
        k0 = harmonic_guard_pick_subharmonic(P, k0, k_lo, harm_max, harm_tol, harm_ratio)
        if k0 <= k_lo:
            return np.nan

    delta = parabolic_peak_interp_logP(P, k0)
    f_est = (k0 + delta) * fs / N
    if not np.isfinite(f_est) or f_est <= 0:
        return np.nan

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


def _col_name(col_idx: int) -> str:
    if col_idx == 1:
        return "total"
    if col_idx == 2:
        return "breath"
    if col_idx == 3:
        return "heart"
    return f"col{col_idx}"


def _pick_signal_by_col(total: float, breath: float, heart: float, col_idx: int) -> float:
    if col_idx == 1:
        return float(total)
    if col_idx == 2:
        return float(breath)
    if col_idx == 3:
        return float(heart)
    return float(heart)


# ===================== STATE =====================
class State:
    def __init__(self):
        self.t0_ms = None
        self.last_t_ms = None

        # shared time
        self.t = deque(maxlen=200000)

        # signals used for FFT (HR and BR)
        self.x_hr = deque(maxlen=200000)
        self.x_br = deque(maxlen=200000)

        # raw channels
        self.t_total  = deque(maxlen=30000)
        self.total    = deque(maxlen=30000)
        self.t_breath = deque(maxlen=30000)
        self.breath   = deque(maxlen=30000)
        self.t_heart  = deque(maxlen=30000)
        self.heart    = deque(maxlen=30000)

        # outputs
        self.hr_t = deque(maxlen=MAXP)
        self.hr   = deque(maxlen=MAXP)
        self.br_t = deque(maxlen=MAXP)
        self.br   = deque(maxlen=MAXP)

        self.fs_t = deque(maxlen=MAXP)
        self.fs   = deque(maxlen=MAXP)

        self.good = 0
        self.bad  = 0
        self.ts_back = 0

        self.lock = threading.Lock()
        self.win_id = 0

        self.csv_path = None
        self.csv_queue = deque(maxlen=CSV_MAX_QUEUE)
        self.csv_lock = threading.Lock()
        self.csv_flush_counter = 0


# ===================== CSV =====================
def _csv_init(st: State):
    if not CSV_ON:
        return
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    st.csv_path = f"ble_fft_hr_br_debug_{ts}.csv"
    with open(st.csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "win_id", "t_center_s", "fs_win_hz",
            "HR_COL", "hr_wi_hz", "hr_wf_hz", "hr_bpm",
            "BR_COL", "br_wi_hz", "br_wf_hz", "br_brpm",
            "t_s", "total", "breath", "heart", "x_used_hr", "x_used_br"
        ])


def _csv_enqueue_window(st: State, win_id, t_center, fs, tw,
                        total_w, breath_w, heart_w, x_hr_w, x_br_w,
                        hr, br):
    if not CSV_ON or st.csv_path is None:
        return

    tw = np.asarray(tw, dtype=float).reshape(-1)
    total_w  = np.asarray(total_w, dtype=float).reshape(-1)
    breath_w = np.asarray(breath_w, dtype=float).reshape(-1)
    heart_w  = np.asarray(heart_w, dtype=float).reshape(-1)
    x_hr_w   = np.asarray(x_hr_w, dtype=float).reshape(-1)
    x_br_w   = np.asarray(x_br_w, dtype=float).reshape(-1)

    n = min(tw.size, total_w.size, breath_w.size, heart_w.size, x_hr_w.size, x_br_w.size)
    if n <= 0:
        return

    rows = []
    for i in range(n):
        rows.append([
            int(win_id), float(t_center), float(fs),
            int(CFG["HR_COL"]), float(CFG["hr_wi"]), float(CFG["hr_wf"]), float(hr) if np.isfinite(hr) else np.nan,
            int(CFG["BR_COL"]), float(CFG["br_wi"]), float(CFG["br_wf"]), float(br) if np.isfinite(br) else np.nan,
            float(tw[i]), float(total_w[i]), float(breath_w[i]), float(heart_w[i]),
            float(x_hr_w[i]), float(x_br_w[i]),
        ])

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


# ===================== PROCESSING =====================
def process_if_possible(st: State, max_windows=2):
    N, hop = int(CFG["N"]), int(CFG["hop"])
    done = 0

    while done < int(max_windows):
        with st.lock:
            n = len(st.t)
            if n < N:
                return

            t_arr   = np.asarray(st.t, dtype=float)
            x_hr_arr = np.asarray(st.x_hr, dtype=float)
            x_br_arr = np.asarray(st.x_br, dtype=float)

            tt_arr  = np.asarray(st.t_total, dtype=float)
            tot_arr = np.asarray(st.total, dtype=float)
            tb_arr  = np.asarray(st.t_breath, dtype=float)
            bre_arr = np.asarray(st.breath, dtype=float)
            th_arr  = np.asarray(st.t_heart, dtype=float)
            hea_arr = np.asarray(st.heart, dtype=float)

        idx_end = n
        idx_start = idx_end - N
        if idx_start < 0:
            return

        tw   = t_arr[idx_start:idx_end]
        x_hr = x_hr_arr[idx_start:idx_end]
        x_br = x_br_arr[idx_start:idx_end]

        fs = _estimate_fs_from_t(tw)
        if not np.isfinite(fs) or fs <= 0:
            return

        t_center = float(tw[0] + 0.5 * (tw[-1] - tw[0]))

        hr = estimate_rate_window(
            x_hr, fs, CFG["hr_f_min_hz"],
            CFG["hr_harm_on"], CFG["hr_harm_max"], CFG["hr_harm_tol"], CFG["hr_harm_ratio"],
            CFG["hr_wi"], CFG["hr_wf"]
        )

        br = estimate_rate_window(
            x_br, fs, CFG["br_f_min_hz"],
            CFG["br_harm_on"], CFG["br_harm_max"], CFG["br_harm_tol"], CFG["br_harm_ratio"],
            CFG["br_wi"], CFG["br_wf"]
        )

        with st.lock:
            st.win_id += 1
            win_id = st.win_id

            st.fs_t.append(t_center); st.fs.append(float(fs))

            if np.isfinite(hr):
                st.hr_t.append(t_center); st.hr.append(float(hr))
            if np.isfinite(br):
                st.br_t.append(t_center); st.br.append(float(br))

        def _window_from_buffers(tref, yref):
            if tref.size < N:
                return None
            i0 = max(0, tref.size - N)
            return tref[i0:], yref[i0:]

        wt = _window_from_buffers(tt_arr, tot_arr)
        wb = _window_from_buffers(tb_arr, bre_arr)
        wh = _window_from_buffers(th_arr, hea_arr)

        if wt and wb and wh:
            tw2, total_w  = wt
            tb2, breath_w = wb
            th2, heart_w  = wh
            _csv_enqueue_window(st, win_id, t_center, fs, tw, total_w, breath_w, heart_w, x_hr, x_br, hr, br)
            _csv_flush(st, force=False)

        # avanço do hop:
        # - se os dois deram NaN, não anda (espera acumular)
        # - se pelo menos um deu valor, anda
        if (not np.isfinite(hr)) and (not np.isfinite(br)):
            return

        with st.lock:
            for _ in range(hop):
                if st.t:
                    st.t.popleft()
                if st.x_hr:
                    st.x_hr.popleft()
                if st.x_br:
                    st.x_br.popleft()

        done += 1


# ===================== GUI =====================
def make_gui():
    plt.ion()
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

    # HR
    l_hr, = ax1.plot([], [], lw=1.5)
    p_hr, = ax1.plot([], [], linestyle="none", marker="o", markersize=4)
    ax1.set_ylabel("HR (bpm)")
    ax1.grid(True)

    # BR
    l_br, = ax2.plot([], [], lw=1.5)
    p_br, = ax2.plot([], [], linestyle="none", marker="o", markersize=4)
    ax2.set_ylabel("BR (brpm)")
    ax2.grid(True)

    # Fs
    l_fs, = ax3.plot([], [], lw=1.0)
    p_fs, = ax3.plot([], [], linestyle="none", marker="o", markersize=4)
    ax3.set_ylabel("Fs_win (Hz)")
    ax3.grid(True)

    # raw
    l_tot, = ax4.plot([], [], lw=1.0)
    l_bre, = ax4.plot([], [], lw=1.0)
    l_hea, = ax4.plot([], [], lw=1.0)
    ax4.set_ylabel("phases (raw)")
    ax4.set_xlabel("t (s)")
    ax4.grid(True)
    ax4.legend(["total", "breath", "heart"], loc="upper right")

    fig.tight_layout()
    return fig, (ax1, ax2, ax3, ax4), (l_hr, p_hr, l_br, p_br, l_fs, p_fs, l_tot, l_bre, l_hea)


def gui_loop_mainthread(st: State, stop_event: threading.Event):
    fig, (ax1, ax2, ax3, ax4), (l_hr, p_hr, l_br, p_br, l_fs, p_fs, l_tot, l_bre, l_hea) = make_gui()
    period = 1.0 / GUI_REFRESH_HZ

    while plt.fignum_exists(fig.number) and (not stop_event.is_set()):
        with st.lock:
            t_hr = np.asarray(st.hr_t, dtype=float)
            y_hr = np.asarray(st.hr, dtype=float)

            t_br = np.asarray(st.br_t, dtype=float)
            y_br = np.asarray(st.br, dtype=float)

            t_fs = np.asarray(st.fs_t, dtype=float)
            y_fs = np.asarray(st.fs, dtype=float)

            tt = np.asarray(st.t_total, dtype=float)
            yt = np.asarray(st.total, dtype=float)

            tb = np.asarray(st.t_breath, dtype=float)
            yb = np.asarray(st.breath, dtype=float)

            th = np.asarray(st.t_heart, dtype=float)
            yh = np.asarray(st.heart, dtype=float)

            good, bad, back = st.good, st.bad, st.ts_back
            csv_path = st.csv_path

        if t_hr.size > 0:
            l_hr.set_data(t_hr, y_hr); p_hr.set_data(t_hr, y_hr)
            ax1.relim(); ax1.autoscale_view()

        if t_br.size > 0:
            l_br.set_data(t_br, y_br); p_br.set_data(t_br, y_br)
            ax2.relim(); ax2.autoscale_view()

        if t_fs.size > 0:
            l_fs.set_data(t_fs, y_fs); p_fs.set_data(t_fs, y_fs)
            ax3.relim(); ax3.autoscale_view()

        if tt.size > 1:
            l_tot.set_data(tt, yt)
        if tb.size > 1:
            l_bre.set_data(tb, yb)
        if th.size > 1:
            l_hea.set_data(th, yh)
        if (tt.size > 1) or (tb.size > 1) or (th.size > 1):
            ax4.relim(); ax4.autoscale_view()

        ax1.set_title(
            f"good={good} bad={bad} ts_back={back} | N={CFG['N']} hop={CFG['hop']} | "
            f"HR_COL={CFG['HR_COL']}({_col_name(CFG['HR_COL'])}) band=[{CFG['hr_wi']:.2f},{CFG['hr_wf']:.2f}]Hz | "
            f"BR_COL={CFG['BR_COL']}({_col_name(CFG['BR_COL'])}) band=[{CFG['br_wi']:.2f},{CFG['br_wf']:.2f}]Hz | "
            f"csv={(csv_path if (CSV_ON and csv_path) else 'OFF')}"
        )

        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(period)

    stop_event.set()
    try:
        plt.close("all")
    except Exception:
        pass


# ===================== THREADS =====================
def proc_loop(st: State, stop_event: threading.Event):
    while not stop_event.is_set():
        try:
            process_if_possible(st, max_windows=2)
        except Exception:
            pass
        time.sleep(0.005)


async def ble_async_main(st: State, stop_event: threading.Event):
    # valida colunas
    if int(CFG.get("HR_COL", 3)) not in (1, 2, 3):
        print(f"[CFG] HR_COL inválido={CFG.get('HR_COL')}. Use 1=total,2=breath,3=heart. Vou usar 3.")
        CFG["HR_COL"] = 3
    if int(CFG.get("BR_COL", 2)) not in (1, 2, 3):
        print(f"[CFG] BR_COL inválido={CFG.get('BR_COL')}. Use 1=total,2=breath,3=heart. Vou usar 2.")
        CFG["BR_COL"] = 2

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
        stop_event.set()
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

        x_hr = _pick_signal_by_col(total, breath, heart, int(CFG.get("HR_COL", 3)))
        x_br = _pick_signal_by_col(total, breath, heart, int(CFG.get("BR_COL", 2)))

        with st.lock:
            st.t.append(t_sec)
            st.x_hr.append(float(x_hr))
            st.x_br.append(float(x_br))

            st.t_total.append(t_sec);  st.total.append(float(total))
            st.t_breath.append(t_sec); st.breath.append(float(breath))
            st.t_heart.append(t_sec);  st.heart.append(float(heart))

            st.good += 1
            good = st.good

        if good % DEBUG_EVERY_N == 0:
            print(
                f"[OK] t_ms={t_ms:9.2f} t={t_sec:8.3f}s | "
                f"total={total:+.3f} breath={breath:+.3f} heart={heart:+.3f} | "
                f"x_hr(col={CFG['HR_COL']}:{_col_name(CFG['HR_COL'])})={x_hr:+.3f} "
                f"x_br(col={CFG['BR_COL']}:{_col_name(CFG['BR_COL'])})={x_br:+.3f}"
            )

    print("Conectando em:", target.address, target.name)
    async with BleakClient(target.address) as client:
        if not client.is_connected:
            print("Falhou conectar.")
            stop_event.set()
            return

        await client.start_notify(CHAR_UUID, on_notify)
        print("Recebendo (feche a janela do plot ou Ctrl+C)...")

        try:
            while not stop_event.is_set():
                await asyncio.sleep(0.2)
        except asyncio.CancelledError:
            pass
        except KeyboardInterrupt:
            stop_event.set()

        try:
            await client.stop_notify(CHAR_UUID)
        except Exception:
            pass

    stop_event.set()


def ble_thread_entry(st: State, stop_event: threading.Event):
    try:
        asyncio.run(ble_async_main(st, stop_event))
    except KeyboardInterrupt:
        stop_event.set()
    except Exception:
        stop_event.set()


# ===================== MAIN =====================
def main():
    st = State()
    _csv_init(st)

    stop_event = threading.Event()

    t_ble = threading.Thread(target=ble_thread_entry, args=(st, stop_event), daemon=True)
    t_ble.start()

    t_proc = threading.Thread(target=proc_loop, args=(st, stop_event), daemon=True)
    t_proc.start()

    try:
        gui_loop_mainthread(st, stop_event)
    except KeyboardInterrupt:
        stop_event.set()

    _csv_flush(st, force=True)


if __name__ == "__main__":
    main()
