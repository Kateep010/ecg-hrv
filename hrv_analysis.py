#!/usr/bin/env python3
"""HRV 分析：讀取 my_ecg_5min.csv（500 Hz），用 neurokit2 做
R 波偵測 → 時域/頻域 HRV → 龐加萊圖，輸出 PNG 與文字報告。

用法:
    python hrv_analysis.py [csv檔] [-o 輸出前綴]
"""

import argparse
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import neurokit2 as nk
import numpy as np
import pandas as pd

FS = 500


def main() -> None:
    ap = argparse.ArgumentParser(description="ECG HRV 分析")
    ap.add_argument("csv", nargs="?", default="my_ecg_5min.csv")
    ap.add_argument("-o", "--prefix", default=None,
                    help="輸出檔名前綴（預設取 CSV 檔名）")
    args = ap.parse_args()
    prefix = args.prefix or args.csv.rsplit(".", 1)[0]

    # ---- 1. 讀檔 + R 波偵測（用原始 adc 欄位） ----
    df = pd.read_csv(args.csv, comment="#")
    raw = df["adc"].astype(float).to_numpy()
    duration = len(raw) / FS

    clean = nk.ecg_clean(raw, sampling_rate=FS)
    _, info = nk.ecg_peaks(clean, sampling_rate=FS)
    rpeaks = info["ECG_R_Peaks"]
    rr = np.diff(rpeaks) / FS * 1000.0  # ms

    # ---- 2. HRV 指標 ----
    hrv_t = nk.hrv_time(rpeaks, sampling_rate=FS)
    hrv_f = nk.hrv_frequency(rpeaks, sampling_rate=FS, normalize=False)
    hrv_n = nk.hrv_nonlinear(rpeaks, sampling_rate=FS)

    rmssd = hrv_t["HRV_RMSSD"].iloc[0]
    sdnn = hrv_t["HRV_SDNN"].iloc[0]
    pnn50 = hrv_t["HRV_pNN50"].iloc[0]
    lf = hrv_f["HRV_LF"].iloc[0]
    hf = hrv_f["HRV_HF"].iloc[0]
    lfhf = hrv_f["HRV_LFHF"].iloc[0]
    sd1 = hrv_n["HRV_SD1"].iloc[0]
    sd2 = hrv_n["HRV_SD2"].iloc[0]

    # ---- 3. 龐加萊圖 + 頻譜 ----
    x, y = rr[:-1], rr[1:]  # RR(n) vs RR(n+1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    ax1.scatter(x, y, s=14, alpha=0.55, color="tab:blue", edgecolors="none")
    lim = [min(x.min(), y.min()) - 30, max(x.max(), y.max()) + 30]
    ax1.plot(lim, lim, "--", color="gray", lw=1, label="RR(n+1)=RR(n)")
    # SD1/SD2 橢圓（沿恆等線方向為 SD2、垂直方向為 SD1）
    m = rr.mean()
    theta = np.linspace(0, 2 * np.pi, 100)
    ex = m + (sd2 * np.cos(theta)) * np.cos(np.pi / 4) - (sd1 * np.sin(theta)) * np.sin(np.pi / 4)
    ey = m + (sd2 * np.cos(theta)) * np.sin(np.pi / 4) + (sd1 * np.sin(theta)) * np.cos(np.pi / 4)
    ax1.plot(ex, ey, color="tab:red", lw=1.8,
             label=f"SD1={sd1:.1f} ms, SD2={sd2:.1f} ms")
    ax1.set_xlim(lim); ax1.set_ylim(lim)
    ax1.set_aspect("equal")
    ax1.set_xlabel("RR(n) (ms)")
    ax1.set_ylabel("RR(n+1) (ms)")
    ax1.set_title("Poincaré plot")
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(alpha=0.3)

    # RR 序列重取樣後做 Welch PSD（與 neurokit2 的頻段一致標示）
    t_rr = rpeaks[1:] / FS
    t_even = np.arange(t_rr[0], t_rr[-1], 0.25)  # 4 Hz
    rr_even = np.interp(t_even, t_rr, rr)
    from scipy.signal import welch
    f, pxx = welch(rr_even - rr_even.mean(), fs=4.0, nperseg=min(1024, len(rr_even)))
    ax2.plot(f, pxx, color="tab:blue", lw=1.2)
    ax2.axvspan(0.04, 0.15, color="tab:orange", alpha=0.18, label="LF 0.04–0.15 Hz")
    ax2.axvspan(0.15, 0.40, color="tab:green", alpha=0.15, label="HF 0.15–0.40 Hz")
    ax2.set_xlim(0, 0.5)
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel("PSD (ms²/Hz)")
    ax2.set_title(f"RR spectrum  (LF/HF = {lfhf:.2f})")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)

    fig.suptitle(f"{args.csv}  |  {duration:.0f} s @ {FS} Hz  |  "
                 f"{len(rpeaks)} beats, mean HR {60000 / rr.mean():.1f} bpm")
    fig.tight_layout()
    png_path = f"{prefix}_hrv.png"
    fig.savefig(png_path, dpi=130)

    # ---- 4. 文字報告 ----
    report = f"""ECG HRV 分析報告
================
產出時間   : {time.strftime('%Y-%m-%d %H:%M:%S')}
資料檔     : {args.csv}
取樣率     : {FS} Hz
記錄長度   : {duration:.1f} s ({len(raw)} 樣本)
電極脫落   : {int(df['lead_off'].sum())} 筆

R 波偵測（neurokit2, 原始 adc 訊號）
------------------------------------
R 波數     : {len(rpeaks)}
RR 間期    : {rr.mean():.1f} ± {rr.std(ddof=1):.1f} ms (範圍 {rr.min():.0f}–{rr.max():.0f})
平均心率   : {60000 / rr.mean():.1f} bpm

時域 HRV
--------
RMSSD      : {rmssd:.1f} ms
SDNN       : {sdnn:.1f} ms
pNN50      : {pnn50:.1f} %

頻域 HRV (Welch)
----------------
LF (0.04–0.15 Hz) : {lf:.1f} ms²
HF (0.15–0.40 Hz) : {hf:.1f} ms²
LF/HF ratio       : {lfhf:.3f}

龐加萊圖（非線性）
------------------
SD1        : {sd1:.1f} ms
SD2        : {sd2:.1f} ms
SD1/SD2    : {sd1 / sd2:.3f}

輸出檔案
--------
圖片       : {png_path}
"""
    txt_path = f"{prefix}_hrv_report.txt"
    with open(txt_path, "w") as fp:
        fp.write(report)
    print(report)
    print(f"已存: {png_path}, {txt_path}")


if __name__ == "__main__":
    main()
