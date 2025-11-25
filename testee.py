import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
#  LEITORES
# ============================================================

def read_phase_csv(path):
    """
    Lê arquivo phases_raw.csv no formato:
        timestamp_ms, total, breath, heart

    Retorna:
        t      : tempo em segundos (0 no início)
        total  : sinal total
        breath : componente respiração
        heart  : componente coração
    """
    data = pd.read_csv(path, comment="#", header=None)

    tpuro  = data.iloc[:, 0].to_numpy()
    total  = data.iloc[:, 1].to_numpy()
    breath = data.iloc[:, 2].to_numpy()
    heart  = data.iloc[:, 3].to_numpy()

    # tempo em segundos começando em zero
    t = (tpuro - tpuro[0]) / 1000.0

    return t, total, breath, heart


def read_radar_csv(path):
    """
    Lê HR do radar no formato:
        right,9,74,,,'1763411195000'
        right,9,70,,,'1763411205032'

    Assumindo:
        col2 = HR (bpm)
        última coluna = timestamp em ms (epoch) entre aspas

    Retorna:
        t  : tempo em segundos (0 no início)
        hr : HR em bpm
    """
    df = pd.read_csv(
        path,
        header=None,
        sep=",",
        quotechar="'",
        engine="python"
    )

    if df.shape[1] < 3:
        raise ValueError(f"Arquivo {path} com colunas insuficientes: {df.shape[1]}")

    hr = df.iloc[:, 2].to_numpy(dtype=float)

    # timestamp é a última coluna
    ts_ms = df.iloc[:, -1].astype("int64").to_numpy()
    t = (ts_ms - ts_ms[0]) / 1000.0

    return t, hr


def read_polar_txt(path):
    """
    Lê arquivo TXT da Polar no formato:
        Phone timestamp;HR [bpm];HRV [ms];Breathing interval [rpm];
        2025-11-17T17:26:16.058;90
        2025-11-17T17:26:17.052;94;12,1
        ...

    Retorna:
        t      : tempo em segundos (0 no início)
        hr     : HR em bpm
        hrv    : HRV em ms (ou None se não existir)
        br_rpm : respiração em rpm (ou None se não existir)
    """
    df = pd.read_csv(
        path,
        sep=";",
        decimal=",",
        engine="python",
        header=0
    )

    time_col = df.columns[0]
    hr_col   = df.columns[1]

    df[time_col] = pd.to_datetime(df[time_col])

    t = (df[time_col] - df[time_col].iloc[0]).dt.total_seconds().to_numpy()
    hr = df[hr_col].to_numpy(dtype=float)

    hrv = None
    br_rpm = None

    if len(df.columns) > 2:
        try:
            hrv = df.iloc[:, 2].to_numpy(dtype=float)
        except Exception:
            hrv = None

    if len(df.columns) > 3:
        try:
            br_rpm = df.iloc[:, 3].to_numpy(dtype=float)
        except Exception:
            br_rpm = None

    return t, hr, hrv, br_rpm


# ============================================================
#  FFT DESLIZANTE (N=128)
# ============================================================

def sliding_fft_band(t, x, N=128, overlap=0.9, band=None):
    """
    FFT deslizante com janela de tamanho N (default 128) e
    banda opcional para procurar o pico.

    t   : vetor de tempo (s)
    x   : sinal
    N   : tamanho da janela (amostras)
    overlap : fração de sobreposição (0.0 a 0.99)
    band: (fmin, fmax) em Hz para restringir a busca do pico

    Retorna:
        t_out : tempo central de cada janela
        f_out : frequência estimada (Hz) em cada janela
        fs    : taxa de amostragem média (Hz)
    """
    x = np.asarray(x)
    t = np.asarray(t)

    dt = np.mean(np.diff(t))
    fs = 1.0 / dt

    hop = max(1, int(N * (1 - overlap)))
    half = N // 2

    f_out = []
    t_out = []

    i = 0
    while i + N <= len(x):
        w = x[i : i + N]
        w = w - np.mean(w)

        W = np.fft.fft(w, n=N)
        Wmag = np.abs(W[:half])

        f = np.fft.fftfreq(N, d=1.0/fs)[:half]

        # tira DC
        if len(Wmag) > 0:
            Wmag[0] = 0.0

        if band is not None:
            fmin, fmax = band
            mask = (f >= fmin) & (f <= fmax)
            if np.any(mask):
                idx_local = np.argmax(Wmag[mask])
                freq = f[mask][idx_local]
            else:
                freq = np.nan
        else:
            idx = np.argmax(Wmag)
            freq = f[idx] if len(f) > 0 else np.nan

        f_out.append(freq)
        # tempo central da janela
        center_idx = i + N // 2
        if center_idx >= len(t):
            center_idx = len(t) - 1
        t_out.append(t[center_idx])

        i += hop

    return np.array(t_out), np.array(f_out), fs


# ============================================================
#  CÁLCULO DE ERROS (RADAR x POLAR)
# ============================================================

def calcular_erros_temporais(t_ref, y_ref, t_est, y_est, nome="HR"):
    """
    Alinha y_est (ex: radar) em t_est com y_ref (ex: Polar) em t_ref via interpolação
    e calcula erros básicos.

    t_ref, y_ref : série de referência (ex: Polar)      -> NÃO interpola
    t_est, y_est : série estimada (ex: Radar via FFT)   -> interpolada em t_ref

    Imprime:
        N pontos usados
        Erro médio (bias)
        MAE
        RMSE
        Desvio padrão do erro
        Correlação
    """
    t_ref = np.asarray(t_ref, dtype=float)
    y_ref = np.asarray(y_ref, dtype=float)
    t_est = np.asarray(t_est, dtype=float)
    y_est = np.asarray(y_est, dtype=float)

    if len(t_ref) == 0 or len(t_est) == 0:
        print(f"[{nome}] Séries vazias, não é possível calcular erro.")
        return

    # usa apenas pontos da referência dentro do intervalo do estimado
    mask_time = (t_ref >= t_est[0]) & (t_ref <= t_est[-1])
    if not np.any(mask_time):
        print(f"[{nome}] Não há sobreposição de tempo entre Radar e Polar.")
        return

    t_sync = t_ref[mask_time]
    y_ref_sync = y_ref[mask_time]

    # interpola estimado no tempo da referência
    y_est_interp = np.interp(t_sync, t_est, y_est)

    # remove NaN (se existirem)
    mask_val = ~np.isnan(y_ref_sync) & ~np.isnan(y_est_interp)
    if np.sum(mask_val) < 2:
        print(f"[{nome}] Dados insuficientes após filtragem para cálculo de erro.")
        return

    y_ref_sync = y_ref_sync[mask_val]
    y_est_interp = y_est_interp[mask_val]

    erro = y_est_interp - y_ref_sync

    mae  = np.mean(np.abs(erro))
    rmse = np.sqrt(np.mean(erro**2))
    bias = np.mean(erro)
    std  = np.std(erro)
    corr = np.corrcoef(y_est_interp, y_ref_sync)[0, 1]

    print(f"\n=== Erros {nome} (Radar vs Polar) ===")
    print(f"N pontos                : {len(erro)}")
    print(f"Erro médio (bias)       : {bias:.3f}")
    print(f"MAE                     : {mae:.3f}")
    print(f"RMSE                    : {rmse:.3f}")
    print(f"Desvio padrão do erro   : {std:.3f}")
    print(f"Correlação              : {corr:.3f}")


# ============================================================
#  MAIN DE EXEMPLO
# ============================================================

if __name__ == "__main__":
    # ==== AJUSTE AQUI OS CAMINHOS DOS ARQUIVOS ====
    phase_path    = r"C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\phases_raw.csv"
    radar_hr_path = r"C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\tes.csv"
    polar_path    = r"C:\Users\eduar\UTFPR\IC\PROJ_ESP32\output\POLARH10.txt"

    # ------------------------------------------------
    # 1) LER PHASE (TOTAL / BREATH / HEART)
    # ------------------------------------------------
    t_phase, total, breath, heart = read_phase_csv(phase_path)

    # ------------------------------------------------
    # 2) FFT DESLIZANTE NO HEART E BREATH (N=128)
    #    bandas:
    #       heart  ~ 0.8 a 3.0 Hz
    #       breath ~ 0.1 a 0.5 Hz
    # ------------------------------------------------
    t_h, f_h, fs_h = sliding_fft_band(
        t_phase, heart,
        N=128,
        overlap=0.9,
        band=(0.8, 3.0)
    )

    t_b, f_b, fs_b = sliding_fft_band(
        t_phase, breath,
        N=128,
        overlap=0.9,
        band=(0.1, 0.5)
    )

    # converter para BPM / RPM
    HR_radar = f_h * 60.0          # bpm
    RR_radar = f_b * 60.0          # rpm (respiração por minuto)

    # ------------------------------------------------
    # 3) (OPCIONAL) LER RADAR HR DIRETO DO CSV TES.CSV
    # ------------------------------------------------
    try:
        t_radar_hr, hr_radar_direct = read_radar_csv(radar_hr_path)
    except Exception as e:
        print("Falha ao ler radar HR CSV:", e)
        t_radar_hr, hr_radar_direct = None, None

    # ------------------------------------------------
    # 4) (OPCIONAL) LER POLAR
    # ------------------------------------------------
    try:
        t_polar, hr_polar, hrv_polar, br_polar = read_polar_txt(polar_path)
    except Exception as e:
        print("Falha ao ler POLAR TXT:", e)
        t_polar, hr_polar, hrv_polar, br_polar = None, None, None, None

    # ------------------------------------------------
    # 5) CÁLCULO DE ERROS (RADAR vs POLAR) – HR
    # ------------------------------------------------
    if t_polar is not None and hr_polar is not None:
        calcular_erros_temporais(t_polar, hr_polar, t_h, HR_radar, nome="HR (bpm)")

    # Se quiser também erro de respiração e tiver br_polar:
    if br_polar is not None:
        calcular_erros_temporais(t_polar, br_polar, t_b, RR_radar, nome="RR (rpm)")

    # ------------------------------------------------
    # 6) PLOTS BÁSICOS
    # ------------------------------------------------

    # HR do radar por phase (FFT)
    plt.figure(figsize=(10, 4))
    plt.scatter(t_h, HR_radar, s=10, label="Radar HR (FFT)")
    plt.title("Radar HR (via phase FFT, N=128)")
    plt.xlabel("Tempo (s)")
    plt.ylabel("HR (bpm)")
    plt.grid(True)

    # RR do radar por phase (FFT)
    plt.figure(figsize=(10, 4))
    plt.scatter(t_b, RR_radar, s=10, color='orange', label="Radar RR (FFT)")
    plt.title("Radar RR (via phase FFT, N=128)")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Respiração (rpm)")
    plt.grid(True)

    # Comparar com POLAR, se disponível
    if t_polar is not None and hr_polar is not None:
        plt.figure(figsize=(10, 4))
        plt.scatter(t_h, HR_radar, s=10, label="Radar HR (FFT)")
        plt.scatter(t_polar, hr_polar, s=10, label="Polar HR", alpha=0.7)
        plt.title("Comparação HR Radar (phase FFT) x Polar")
        plt.xlabel("Tempo (s)")
        plt.ylabel("HR (bpm)")
        plt.legend()
        plt.grid(True)

    # HR direto do CSV do radar, se disponível
    if t_radar_hr is not None and hr_radar_direct is not None:
        plt.figure(figsize=(10, 4))
        plt.scatter(t_radar_hr, hr_radar_direct, s=10)
        plt.title("Radar HR direto do CSV (tes.csv)")
        plt.xlabel("Tempo (s)")
        plt.ylabel("HR (bpm)")
        plt.grid(True)

    plt.show()
