import asyncio
import struct
import time
import threading
from collections import deque

import numpy as np
import matplotlib.pyplot as plt
from bleak import BleakScanner, BleakClient


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

GUI_REFRESH_HZ  = 10.0
MAXP            = 1500
DEBUG_EVERY_N   = 50


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

        self.fft_t = deque(maxlen=MAXP)
        self.fft_ms = deque(maxlen=MAXP)

        self.good = 0
        self.bad  = 0
        self.ts_back = 0

        self.lock = threading.Lock()


def process_if_possible(st: State):
    N, hop = CFG["N"], CFG["hop"]

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
    xw = x_arr[idx_start:idx_end]

    dt = np.diff(tw)
    ok = np.isfinite(dt) & (dt > 0)
    if np.count_nonzero(ok) < max(3, N//4):
        return
    fs = 1.0 / np.mean(dt[ok])

    t_center = float(tw[0] + 0.5*(tw[-1] - tw[0]))

    t0 = time.perf_counter()
    hr = estimate_hr_window(
        xw, fs, CFG["f_min_hz"],
        CFG["harm_on"], CFG["harm_max"], CFG["harm_tol"], CFG["harm_ratio"]
    )
    t1 = time.perf_counter()
    dt_ms = 1000.0 * (t1 - t0)

    if not np.isfinite(hr):
        return

    with st.lock:
        st.hr_t.append(t_center); st.hr.append(float(hr))
        st.fs_t.append(t_center); st.fs.append(float(fs))
        st.fft_t.append(t_center); st.fft_ms.append(float(dt_ms))

        for _ in range(hop):
            if st.x:
                st.x.popleft()
                st.t.popleft()


def make_gui():
    plt.ion()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7), sharex=True)

    l_hr, = ax1.plot([], [], lw=1.5)
    ax1.set_ylabel("HR (bpm)")
    ax1.grid(True)

    l_fs, = ax2.plot([], [], lw=1.0)
    ax2.set_ylabel("Fs_win (Hz)")
    ax2.grid(True)

    l_fft, = ax3.plot([], [], lw=1.0)
    ax3.set_ylabel("FFT time (ms)")
    ax3.set_xlabel("t (s)")
    ax3.grid(True)

    fig.tight_layout()
    return fig, (ax1, ax2, ax3), (l_hr, l_fs, l_fft)


def gui_loop(st: State):
    fig, (ax1, ax2, ax3), (l_hr, l_fs, l_fft) = make_gui()
    period = 1.0 / GUI_REFRESH_HZ

    while plt.fignum_exists(fig.number):
        with st.lock:
            t_hr = np.asarray(st.hr_t, dtype=float)
            y_hr = np.asarray(st.hr, dtype=float)

            t_fs = np.asarray(st.fs_t, dtype=float)
            y_fs = np.asarray(st.fs, dtype=float)

            t_ft = np.asarray(st.fft_t, dtype=float)
            y_ft = np.asarray(st.fft_ms, dtype=float)

            good, bad, back = st.good, st.bad, st.ts_back

        if t_hr.size > 0:
            l_hr.set_data(t_hr, y_hr)
            ax1.relim(); ax1.autoscale_view()

        if t_fs.size > 0:
            l_fs.set_data(t_fs, y_fs)
            ax2.relim(); ax2.autoscale_view()

        if t_ft.size > 0:
            l_fft.set_data(t_ft, y_ft)
            ax3.relim(); ax3.autoscale_view()

        ax1.set_title(f"good={good} bad={bad} ts_back={back} | N={CFG['N']} hop={CFG['hop']}")

        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(period)


async def main():
    st = State()

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

    gui_thread = threading.Thread(target=gui_loop, args=(st,), daemon=True)
    gui_thread.start()

    def on_notify(sender: int, data: bytearray):
        nonlocal st

        rec, err = parse_pkt_19(data)
        if rec is None:
            st.bad += 1
            if st.bad % 20 == 0:
                print(f"[BAD] {err} hex={data.hex()}")
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

        process_if_possible(st)

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
        print("Parou.")


if __name__ == "__main__":
    asyncio.run(main())
