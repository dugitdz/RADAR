import asyncio
import struct
import time
import threading
from collections import deque
from queue import Queue, Empty
import os
import csv
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from bleak import BleakScanner, BleakClient

from scipy.signal import cheby2, sosfilt, sosfilt_zi


# ===================== BLE CONFIG =====================
DEVICE_NAME = "ESP32C3_IC"
CHAR_UUID   = "00001234-0000-1000-8000-00805f9b34fb"

# ===================== CSV CONFIG =====================
CSV_DIR  = r"C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\HR_BR"
CSV_ON   = 1
CSV_FLUSH_EVERY = 50  # flush a cada N linhas

# ===================== CONFIG =====================
CFG = dict(
    # quais canais do pacote usar (1=total,2=breath,3=heart)
    HR_COL=1,
    BR_COL=2,

    # tracker (N/hop separados)
    N=256,        # HR
    hop=32,       # HR
    N_BR=512,     # BR
    hop_BR=64,    # BR

    # "fs" assumido pelo tracker (sem resample)
    FS_TARGET=50.0,

    # segmentação por gap (segundos)
    GAP_THR=0.6,

    # --- HR bandpass causal (cheby2, ordem PAR) ---
    HR_FILT_ORD=4,
    HR_FILT_WI=0.595,
    HR_FILT_WF=3.05,
    HR_FILT_RS=30.50,

    # --- BR bandpass causal (cheby2, ordem PAR) ---
    BR_FILT_ORD=4,
    BR_FILT_WI=0.10,
    BR_FILT_WF=0.50,
    BR_FILT_RS=30.0,

    # prior HR (Hz)
    HR_PRIOR_INIT_BPM=80.0,
    HR_PRIOR_MIN_BPM=35.0,
    HR_PRIOR_MAX_BPM=240.0,
    HR_PRIOR_SIGMA_HZ=0.40,
    HR_PRIOR_ALPHA=0.05,
    HR_EMA_TAU_SEC=4.0,

    # prior BR (Hz)
    BR_PRIOR_INIT_BRPM=15.0,
    BR_PRIOR_MIN_BRPM=4.0,
    BR_PRIOR_MAX_BRPM=40.0,
    BR_PRIOR_SIGMA_HZ=0.10,
    BR_PRIOR_ALPHA=0.10,
    BR_EMA_TAU_SEC=6.0,

    # FFT min freq (evita bins muito baixos)
    F_MIN_HZ=0.05,
)

GUI_REFRESH_HZ = 20.0
MAXP           = 2000
DEBUG_EVERY_N  = 80


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
    return (float(t_ms), float(total), float(breath), float(heart)), None

def _pick_signal_by_col(total: float, breath: float, heart: float, col_idx: int) -> float:
    if col_idx == 1:
        return float(total)
    if col_idx == 2:
        return float(breath)
    return float(heart)

def _col_name(i: int) -> str:
    return {1: "total", 2: "resp", 3: "heart"}.get(int(i), f"col{i}")


# ===================== DSP (realtime-friendly) =====================
def design_sos_cheby2(fs, order, wi, wf, rs):
    fs = float(fs)
    order = int(order)
    wi = float(wi); wf = float(wf); rs = float(rs)

    if (not np.isfinite(fs)) or fs <= 0:
        return None, None, "bad_fs"
    if (order % 2) != 0 or order < 2:
        return None, None, "order_must_be_even"
    if (not np.isfinite(wi)) or (not np.isfinite(wf)) or wi <= 0 or wf <= wi:
        return None, None, "bad_wi_wf"
    if wf >= fs/2:
        return None, None, "wf>=nyq"

    try:
        sos = cheby2(order // 2, rs, [wi/(fs/2), wf/(fs/2)], btype="bandpass", output="sos")
        zi  = sosfilt_zi(sos)
        return sos, zi, None
    except Exception as e:
        return None, None, f"design_fail:{type(e).__name__}"

def parabolic_peak(P, k):
    if k <= 0 or k >= (len(P) - 1):
        return 0.0
    val = np.log(P[k-1:k+2] + 1e-20)
    den = val[0] - 2*val[1] + val[2]
    if den == 0 or (not np.isfinite(den)):
        return 0.0
    return float(np.clip(0.5 * (val[0] - val[2]) / den, -0.75, 0.75))


class RealtimeTracker:
    def __init__(self, fs, N, hop, f_min_hz,
                 prior_init_per_min, prior_min_per_min, prior_max_per_min,
                 prior_sigma_hz, prior_alpha, ema_tau_sec):

        self.fs = float(fs)
        self.N = int(N)
        self.hop = int(hop)
        self.f_min_hz = float(f_min_hz)

        # janela flattop (periodic)
        n = np.arange(self.N, dtype=float)
        self.win = (0.21557895
                    - 0.41663158*np.cos(2*np.pi*n/self.N)
                    + 0.277263158*np.cos(4*np.pi*n/self.N)
                    - 0.083578947*np.cos(6*np.pi*n/self.N)
                    + 0.006947368*np.cos(8*np.pi*n/self.N))

        self.k_min = max(1, int(np.ceil(self.f_min_hz * self.N / self.fs)) - 1)
        self.f_bins = np.fft.rfftfreq(self.N, d=1.0/self.fs)

        self.prior_min = float(prior_min_per_min) / 60.0
        self.prior_max = float(prior_max_per_min) / 60.0
        self.f_prev = np.clip(float(prior_init_per_min) / 60.0, self.prior_min, self.prior_max)

        self.prior_sigma_hz = float(prior_sigma_hz)
        self.prior_alpha = float(prior_alpha)

        if ema_tau_sec and float(ema_tau_sec) > 0:
            self.ema_alpha = 1.0 - np.exp(-(self.hop/self.fs) / float(ema_tau_sec))
        else:
            self.ema_alpha = 1.0
        self.ema_state = None

        self.buf = deque(maxlen=self.N)
        self.tbuf = deque(maxlen=self.N)
        self._cnt = 0

        self.t_out = deque(maxlen=MAXP)
        self.y_out = deque(maxlen=MAXP)

    def reset(self):
        self.buf.clear()
        self.tbuf.clear()
        self._cnt = 0
        self.ema_state = None

    def push(self, t, x):
        self.buf.append(float(x))
        self.tbuf.append(float(t))
        self._cnt += 1

        if len(self.buf) < self.N:
            return None

        if self._cnt < self.hop:
            return None
        self._cnt = 0

        xw = np.array(self.buf, dtype=float) * self.win
        xw = xw - np.mean(xw)

        X = np.fft.rfft(xw, n=self.N)
        P = np.abs(X)**2

        k_srch_min = max(self.k_min, int(np.ceil(self.prior_min * self.N / self.fs)))
        k_srch_max = min(len(P) - 1, int(np.floor(self.prior_max * self.N / self.fs)))

        if k_srch_max > k_srch_min and np.isfinite(self.prior_sigma_hz) and self.prior_sigma_hz > 0:
            W = np.exp(-0.5 * ((self.f_bins - self.f_prev) / self.prior_sigma_hz)**2)
            k0 = k_srch_min + int(np.argmax((P * W)[k_srch_min:k_srch_max + 1]))
        else:
            k0 = self.k_min + int(np.argmax(P[self.k_min:]))

        f_est = (k0 + parabolic_peak(P, k0)) * self.fs / self.N
        if (not np.isfinite(f_est)) or f_est <= 0:
            return None

        y_inst_per_min = 60.0 * f_est
        if self.ema_state is None:
            self.ema_state = y_inst_per_min
        else:
            self.ema_state = (1 - self.ema_alpha) * self.ema_state + self.ema_alpha * y_inst_per_min

        f_clipped = np.clip(f_est, self.prior_min, self.prior_max)
        self.f_prev = (1 - self.prior_alpha) * self.f_prev + self.prior_alpha * f_clipped

        t_center = self.tbuf[0] + 0.5 * (self.N - 1) / self.fs
        self.t_out.append(float(t_center))
        self.y_out.append(float(self.ema_state))
        return float(self.ema_state)


# ===================== CSV WRITER =====================
def csv_writer_thread(st, q: Queue, out_path: str):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    n_written = 0

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=",")
        # timestamp,total,resp,heart,BR,HR
        w.writerow(["timestamp", "total", "resp", "heart", "BR", "HR"])

        while (not st.stop_event.is_set()) or (not q.empty()):
            try:
                row = q.get(timeout=0.2)
            except Empty:
                continue

            w.writerow(row)
            n_written += 1
            if (CSV_FLUSH_EVERY > 0) and (n_written % CSV_FLUSH_EVERY == 0):
                f.flush()

        f.flush()


# ===================== STATE =====================
class State:
    def __init__(self):
        self.lock = threading.Lock()
        self.stop_event = threading.Event()

        self.t0_ms = None
        self.last_t_sec = None

        self.good = 0
        self.bad = 0
        self.gaps = 0
        self.ts_back = 0

        # plot buffers
        self.t_raw    = deque(maxlen=30000)
        self.heart_r  = deque(maxlen=30000)
        self.breath_r = deque(maxlen=30000)

        self.heart_f  = deque(maxlen=30000)
        self.breath_f = deque(maxlen=30000)

        self.t_hr = deque(maxlen=MAXP)
        self.hr   = deque(maxlen=MAXP)

        self.t_br = deque(maxlen=MAXP)
        self.br   = deque(maxlen=MAXP)

        # últimas estimativas (para logar em toda linha do CSV)
        self.last_hr = float("nan")
        self.last_br = float("nan")

        self._init_dsp_objects()

    def _init_dsp_objects(self):
        fs = CFG["FS_TARGET"]

        sos_h, zi_h, err_h = design_sos_cheby2(
            fs, CFG["HR_FILT_ORD"], CFG["HR_FILT_WI"], CFG["HR_FILT_WF"], CFG["HR_FILT_RS"]
        )
        if sos_h is None:
            raise SystemExit(f"HR filter error: {err_h}")

        sos_b, zi_b, err_b = design_sos_cheby2(
            fs, CFG["BR_FILT_ORD"], CFG["BR_FILT_WI"], CFG["BR_FILT_WF"], CFG["BR_FILT_RS"]
        )
        if sos_b is None:
            raise SystemExit(f"BR filter error: {err_b}")

        self.sos_hr = sos_h
        self.zi_hr0 = zi_h.copy()
        self.zi_hr  = zi_h.copy()

        self.sos_br = sos_b
        self.zi_br0 = zi_b.copy()
        self.zi_br  = zi_b.copy()

        fmin = CFG["F_MIN_HZ"]

        self.trk_hr = RealtimeTracker(
            fs=fs, N=CFG["N"], hop=CFG["hop"], f_min_hz=fmin,
            prior_init_per_min=CFG["HR_PRIOR_INIT_BPM"],
            prior_min_per_min=CFG["HR_PRIOR_MIN_BPM"],
            prior_max_per_min=CFG["HR_PRIOR_MAX_BPM"],
            prior_sigma_hz=CFG["HR_PRIOR_SIGMA_HZ"],
            prior_alpha=CFG["HR_PRIOR_ALPHA"],
            ema_tau_sec=CFG["HR_EMA_TAU_SEC"],
        )

        self.trk_br = RealtimeTracker(
            fs=fs, N=CFG["N_BR"], hop=CFG["hop_BR"], f_min_hz=fmin,
            prior_init_per_min=CFG["BR_PRIOR_INIT_BRPM"],
            prior_min_per_min=CFG["BR_PRIOR_MIN_BRPM"],
            prior_max_per_min=CFG["BR_PRIOR_MAX_BRPM"],
            prior_sigma_hz=CFG["BR_PRIOR_SIGMA_HZ"],
            prior_alpha=CFG["BR_PRIOR_ALPHA"],
            ema_tau_sec=CFG["BR_EMA_TAU_SEC"],
        )

    def reset_segment(self, x_hr_init=0.0, x_br_init=0.0):
        self.zi_hr = self.zi_hr0 * float(x_hr_init)
        self.zi_br = self.zi_br0 * float(x_br_init)
        self.trk_hr.reset()
        self.trk_br.reset()
        # mantém last_hr/last_br (você pode zerar aqui se quiser)
        # self.last_hr = float("nan")
        # self.last_br = float("nan")


# ===================== GUI =====================
def make_gui():
    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 9), sharex=True)

    l_hr, = ax1.plot([], [], lw=1.5)
    p_hr, = ax1.plot([], [], linestyle="none", marker="o", markersize=4)
    ax1.set_ylabel("HR (bpm)")
    ax1.grid(True)

    l_br, = ax2.plot([], [], lw=1.5)
    p_br, = ax2.plot([], [], linestyle="none", marker="o", markersize=4)
    ax2.set_ylabel("BR (brpm)")
    ax2.grid(True)

    l_h_r, = ax3.plot([], [], lw=1.0)
    l_b_r, = ax3.plot([], [], lw=1.0)
    l_h_f, = ax3.plot([], [], lw=1.5)
    l_b_f, = ax3.plot([], [], lw=1.5)

    ax3.set_ylabel("signals")
    ax3.set_xlabel("t (s)")
    ax3.grid(True)
    ax3.legend(["heart_raw", "resp_raw", "heart_filt", "resp_filt"], loc="upper right")

    fig.tight_layout()
    return fig, (ax1, ax2, ax3), (l_hr, p_hr, l_br, p_br, l_h_r, l_b_r, l_h_f, l_b_f)


def gui_loop_mainthread(st: State):
    fig, (ax1, ax2, ax3), (l_hr, p_hr, l_br, p_br, l_h_r, l_b_r, l_h_f, l_b_f) = make_gui()
    period = 1.0 / float(GUI_REFRESH_HZ)

    while plt.fignum_exists(fig.number) and (not st.stop_event.is_set()):
        with st.lock:
            t_hr = np.asarray(st.t_hr, dtype=float)
            y_hr = np.asarray(st.hr, dtype=float)

            t_br = np.asarray(st.t_br, dtype=float)
            y_br = np.asarray(st.br, dtype=float)

            t   = np.asarray(st.t_raw, dtype=float)
            h_r = np.asarray(st.heart_r, dtype=float)
            b_r = np.asarray(st.breath_r, dtype=float)
            h_f = np.asarray(st.heart_f, dtype=float)
            b_f = np.asarray(st.breath_f, dtype=float)

            good, bad, gaps, back = st.good, st.bad, st.gaps, st.ts_back

        if t_hr.size > 0:
            l_hr.set_data(t_hr, y_hr); p_hr.set_data(t_hr, y_hr)
            ax1.relim(); ax1.autoscale_view()

        if t_br.size > 0:
            l_br.set_data(t_br, y_br); p_br.set_data(t_br, y_br)
            ax2.relim(); ax2.autoscale_view()

        if t.size > 1:
            l_h_r.set_data(t, h_r)
            l_b_r.set_data(t, b_r)
            l_h_f.set_data(t, h_f)
            l_b_f.set_data(t, b_f)
            ax3.relim(); ax3.autoscale_view()

        ax1.set_title(
            f"good={good} bad={bad} gaps={gaps} ts_back={back} | "
            f"HR: N={CFG['N']} hop={CFG['hop']} | BR: N_BR={CFG['N_BR']} hop_BR={CFG['hop_BR']} | fs={CFG['FS_TARGET']:.2f} | "
            f"HR_COL={CFG['HR_COL']}({_col_name(CFG['HR_COL'])}) BP=[{CFG['HR_FILT_WI']:.3f},{CFG['HR_FILT_WF']:.3f}]Hz | "
            f"BR_COL={CFG['BR_COL']}({_col_name(CFG['BR_COL'])}) BP=[{CFG['BR_FILT_WI']:.3f},{CFG['BR_FILT_WF']:.3f}]Hz"
        )

        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(period)

    st.stop_event.set()
    try:
        plt.close("all")
    except Exception:
        pass


# ===================== BLE THREAD (async) =====================
async def ble_async_main(st: State, csv_q: Queue):
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
        st.stop_event.set()
        return

    def on_notify(sender: int, data: bytearray):
        rec, err = parse_pkt_19(data)
        if rec is None:
            with st.lock:
                st.bad += 1
            return

        t_ms, total, resp, heart = rec

        with st.lock:
            if st.t0_ms is None:
                st.t0_ms = float(t_ms)

        t_sec = (float(t_ms) - float(st.t0_ms)) / 1000.0

        with st.lock:
            if st.last_t_sec is not None and t_sec < st.last_t_sec:
                st.ts_back += 1

        gap = False
        with st.lock:
            if st.last_t_sec is not None and (t_sec - st.last_t_sec) > float(CFG["GAP_THR"]):
                st.gaps += 1
                gap = True
            st.last_t_sec = t_sec

        x_hr_raw = _pick_signal_by_col(total, resp, heart, int(CFG["HR_COL"]))
        x_br_raw = _pick_signal_by_col(total, resp, heart, int(CFG["BR_COL"]))

        if gap:
            with st.lock:
                st.reset_segment(x_hr_init=x_hr_raw, x_br_init=x_br_raw)

        with st.lock:
            y_hr, st.zi_hr = sosfilt(st.sos_hr, [x_hr_raw], zi=st.zi_hr)
            y_br, st.zi_br = sosfilt(st.sos_br, [x_br_raw], zi=st.zi_br)
            y_hr = float(y_hr[0])
            y_br = float(y_br[0])

            hr_est = st.trk_hr.push(t_sec, y_hr)
            br_est = st.trk_br.push(t_sec, y_br)

            # buffers de plot
            st.t_raw.append(t_sec)
            st.heart_r.append(float(x_hr_raw))
            st.breath_r.append(float(x_br_raw))
            st.heart_f.append(y_hr)
            st.breath_f.append(y_br)

            if hr_est is not None:
                st.t_hr.append(st.trk_hr.t_out[-1])
                st.hr.append(st.trk_hr.y_out[-1])
                st.last_hr = float(st.trk_hr.y_out[-1])

            if br_est is not None:
                st.t_br.append(st.trk_br.t_out[-1])
                st.br.append(st.trk_br.y_out[-1])
                st.last_br = float(st.trk_br.y_out[-1])

            # >>> CSV: 1 linha por pacote BLE <<<
            if CSV_ON:
                csv_q.put((
                    float(t_ms),           # timestamp (ms) como no pacote
                    float(total),
                    float(resp),
                    float(heart),
                    float(st.last_br),     # última BR conhecida (NaN se ainda não existe)
                    float(st.last_hr),     # última HR conhecida (NaN se ainda não existe)
                ))

            st.good += 1
            good = st.good

        if good % DEBUG_EVERY_N == 0:
            with st.lock:
                _hr = st.last_hr
                _br = st.last_br
            print(f"[OK] t={t_sec:8.3f}s | HR={_hr:6.2f} | BR={_br:6.2f}")

    print("Conectando em:", target.address, target.name)
    async with BleakClient(target.address) as client:
        if not client.is_connected:
            print("Falhou conectar.")
            st.stop_event.set()
            return

        await client.start_notify(CHAR_UUID, on_notify)
        print("Recebendo (feche a janela do plot ou Ctrl+C)...")

        try:
            while not st.stop_event.is_set():
                await asyncio.sleep(0.2)
        except KeyboardInterrupt:
            st.stop_event.set()
        finally:
            try:
                await client.stop_notify(CHAR_UUID)
            except Exception:
                pass


def ble_thread_entry(st: State, csv_q: Queue):
    try:
        asyncio.run(ble_async_main(st, csv_q))
    except KeyboardInterrupt:
        st.stop_event.set()
    except Exception:
        st.stop_event.set()


# ===================== MAIN =====================
def main():
    # valida cfg básica
    if int(CFG["N"]) < 32 or int(CFG["hop"]) < 1 or int(CFG["hop"]) > int(CFG["N"]):
        raise SystemExit("CFG inválido (HR): verifique N/hop")
    if int(CFG["N_BR"]) < 32 or int(CFG["hop_BR"]) < 1 or int(CFG["hop_BR"]) > int(CFG["N_BR"]):
        raise SystemExit("CFG inválido (BR): verifique N_BR/hop_BR")
    if int(CFG["HR_COL"]) not in (1, 2, 3) or int(CFG["BR_COL"]) not in (1, 2, 3):
        raise SystemExit("CFG inválido: HR_COL/BR_COL devem ser 1,2,3")

    st = State()

    # CSV (writer thread)
    csv_q = Queue(maxsize=20000)
    if CSV_ON:
        os.makedirs(CSV_DIR, exist_ok=True)
        fname = f"hr_br_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        csv_path = os.path.join(CSV_DIR, fname)
        t_csv = threading.Thread(target=csv_writer_thread, args=(st, csv_q, csv_path), daemon=True)
        t_csv.start()
        print("CSV:", csv_path)

    # BLE thread
    t_ble = threading.Thread(target=ble_thread_entry, args=(st, csv_q), daemon=True)
    t_ble.start()

    # GUI no main thread
    try:
        gui_loop_mainthread(st)
    except KeyboardInterrupt:
        st.stop_event.set()


if __name__ == "__main__":
    main()