#!/usr/bin/env python3
"""產生單頁 HRV 完整報告：ECG+R 波、RR 間期趨勢、時域指標、頻域 LF/HF。

用法:
    python make_hrv_report.py my_ecg_5min.csv
輸出:
    <prefix>_full_report.png / <prefix>_full_report.txt
"""

import argparse
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import neurokit2 as nk
import numpy as np
import pandas as pd
from scipy.signal import welch

FS = 500


def main() -> None:
    ap = argparse.ArgumentParser(description="HRV 完整報告")
    ap.add_argument("csv")
    ap.add_argument("-o", "--prefix", default=None)
    args = ap.parse_args()
    prefix = args.prefix or args.csv.rsplit(".", 1)[0]

    # ---- 分析 ----
    df = pd.read_csv(args.csv, comment="#")
    raw = df["adc"].astype(float).to_numpy()
    duration = len(raw) / FS
    clean = nk.ecg_clean(raw, sampling_rate=FS)
    _, info = nk.ecg_peaks(clean, sampling_rate=FS, correct_artifacts=True)
    rp = np.array(info["ECG_R_Peaks"])
    rr = np.diff(rp) / FS * 1000.0
    t_rr = rp[1:] / FS

    ht = nk.hrv_time(rp, sampling_rate=FS)
    hf_ = nk.hrv_frequency(rp, sampling_rate=FS, normalize=False)
    rmssd = ht["HRV_RMSSD"].iloc[0]
    sdnn = ht["HRV_SDNN"].iloc[0]
    pnn50 = ht["HRV_pNN50"].iloc[0]
    lf = hf_["HRV_LF"].iloc[0]
    hf = hf_["HRV_HF"].iloc[0]
    lfhf = hf_["HRV_LFHF"].iloc[0]
    mean_hr = 60000 / rr.mean()

    # ---- 圖：2x2 版面 ----
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1, 1])

    # (1) ECG 全程 + R 波（上排跨兩欄）
    ax1 = fig.add_subplot(gs[0, :])
    t = np.arange(len(clean)) / FS
    ax1.plot(t, clean, lw=0.4, color="tab:blue")
    ax1.plot(rp / FS, clean[rp], "rv", ms=4, label=f"R peaks (n={len(rp)})")
    ax1.set_xlim(0, duration)
    ax1.set_ylabel("ECG (cleaned, a.u.)")
    ax1.set_title("ECG with R-peak detection")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(alpha=0.25)

    # (2) RR 間期趨勢（tachogram + 30 拍滾動中位數）
    ax2 = fig.add_subplot(gs[1, :], sharex=ax1)
    ax2.plot(t_rr, rr, "o-", ms=2.5, lw=0.7, color="tab:red", alpha=0.7,
             label="RR(n)")
    roll = pd.Series(rr).rolling(30, center=True, min_periods=5).median()
    ax2.plot(t_rr, roll, lw=2, color="darkred", label="30-beat rolling median")
    ax2.set_ylabel("RR interval (ms)")
    ax2.set_xlabel("Time (s)")
    ax2.set_title(f"RR interval trend  (mean {rr.mean():.0f} ms ≈ {mean_hr:.1f} bpm)")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(alpha=0.25)

    # (3) 頻域 PSD + LF/HF
    ax3 = fig.add_subplot(gs[2, 0])
    t_even = np.arange(t_rr[0], t_rr[-1], 0.25)
    rr_even = np.interp(t_even, t_rr, rr)
    f, pxx = welch(rr_even - rr_even.mean(), fs=4.0,
                   nperseg=min(1024, len(rr_even)))
    ax3.plot(f, pxx, lw=1.2, color="tab:blue")
    ax3.fill_between(f, pxx, where=(f >= 0.04) & (f < 0.15),
                     color="tab:orange", alpha=0.4,
                     label=f"LF = {lf:.0f} ms²")
    ax3.fill_between(f, pxx, where=(f >= 0.15) & (f <= 0.40),
                     color="tab:green", alpha=0.4,
                     label=f"HF = {hf:.0f} ms²")
    ax3.set_xlim(0, 0.5)
    ax3.set_xlabel("Frequency (Hz)")
    ax3.set_ylabel("PSD (ms²/Hz)")
    ax3.set_title(f"Frequency domain  (LF/HF = {lfhf:.2f})")
    ax3.legend(fontsize=9)
    ax3.grid(alpha=0.25)

    # (4) 指標摘要
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.axis("off")
    summary = (
        f"Recording   {args.csv}\n"
        f"Duration    {duration:.1f} s @ {FS} Hz\n"
        f"Beats       {len(rp)}  (lead-off: {int(df['lead_off'].sum())})\n"
        f"Mean HR     {mean_hr:.1f} bpm\n"
        f"\n─ Time domain ─────────────\n"
        f"RMSSD       {rmssd:.1f} ms\n"
        f"SDNN        {sdnn:.1f} ms\n"
        f"pNN50       {pnn50:.1f} %\n"
        f"\n─ Frequency domain ────────\n"
        f"LF          {lf:.1f} ms²\n"
        f"HF          {hf:.1f} ms²\n"
        f"LF/HF       {lfhf:.3f}"
    )
    ax4.text(0.05, 0.95, summary, transform=ax4.transAxes, va="top",
             family="monospace", fontsize=11)

    fig.suptitle(f"HRV Report — {args.csv}   ({time.strftime('%Y-%m-%d %H:%M')})",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    png = f"{prefix}_full_report.png"
    fig.savefig(png, dpi=130)

    # ---- 文字報告 ----
    txt = f"{prefix}_full_report.txt"
    with open(txt, "w") as fp:
        fp.write(f"""HRV 完整報告
============
產出時間 : {time.strftime('%Y-%m-%d %H:%M:%S')}
資料檔   : {args.csv}
長度     : {duration:.1f} s @ {FS} Hz
R 波     : {len(rp)} 個（含 artifact 校正）
電極脫落 : {int(df['lead_off'].sum())} 筆
平均心率 : {mean_hr:.1f} bpm
RR 間期  : {rr.mean():.1f} ± {rr.std(ddof=1):.1f} ms（範圍 {rr.min():.0f}–{rr.max():.0f}）

時域指標
--------
RMSSD  : {rmssd:.1f} ms
SDNN   : {sdnn:.1f} ms
pNN50  : {pnn50:.1f} %

頻域指標 (Welch)
----------------
LF (0.04–0.15 Hz) : {lf:.1f} ms²
HF (0.15–0.40 Hz) : {hf:.1f} ms²
LF/HF             : {lfhf:.3f}

輸出: {png}
""")
    print(f"完成: {png}, {txt}")


if __name__ == "__main__":
    main()
