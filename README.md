# ecg-hrv — ECG 心率變異性（HRV）分析

用 Arduino UNO R4 WiFi + AD8232 心電模組自錄 ECG，以 NeuroKit2 做 R 波偵測與 HRV 分析。

## 硬體與資料來源

| 項目 | 規格 |
|---|---|
| 開發板 | Arduino UNO R4 WiFi（`arduino:renesas_uno:unor4wifi`） |
| 感測器 | AD8232 單導程 ECG 模組（OUTPUT→A0、LO+→D10、LO−→D11、3.3V 供電） |
| 取樣率 | 500 Hz，14-bit ADC（0–16383） |
| 韌體/收錄工具 | 見姊妹專案 `my-ecg`（`ecg.ino` 韌體與 `ecg_live.py` 即時顯示/錄製） |

## 資料格式

CSV 欄位：`sample, adc, filtered, lead_off`

- `adc`：原始 ADC 值（分析都用這一欄）
- `filtered`：錄製端即時濾波值（0.5 Hz 高通 → 60 Hz 陷波 → 40 Hz 低通）
- `lead_off`：電極脫落旗標（0/1）

| 資料檔 | 錄製日期 | 長度 |
|---|---|---|
| `my_ecg_5min.csv` | 2026-08-19 | 333 s |
| `my_ecg_5min_0820.csv` | 2026-08-20 | 340 s |

## 使用方式

```bash
conda activate bioagent   # Python 3.11 + neurokit2, numpy, matplotlib, scipy, pandas

# 基本分析：時域/頻域 HRV + 龐加萊圖 + 頻譜
python hrv_analysis.py my_ecg_5min.csv

# 完整單頁報告：ECG+R 波、RR 趨勢、頻譜 LF/HF、龐加萊圖、指標摘要
python make_hrv_report.py my_ecg_5min.csv
```

輸出為 `<資料檔名>_hrv.png/.txt` 與 `<資料檔名>_full_report.png/.txt`。

## 分析流程

1. `neurokit2.ecg_clean` 清理原始 `adc` 訊號
2. `neurokit2.ecg_peaks` R 波偵測（報告版開啟 `correct_artifacts=True`）
3. 時域指標：RMSSD、SDNN、pNN50
4. 頻域指標：RR 序列 4 Hz 等距重取樣 → Welch PSD → LF (0.04–0.15 Hz)、HF (0.15–0.40 Hz)、LF/HF
5. 龐加萊圖：RR(n) vs RR(n+1)，SD1/SD2 橢圓

## 結果摘要

| 指標 | 8/19 | 8/20 |
|---|---|---|
| 平均心率 | 66.1 bpm | 68.1 bpm |
| RMSSD | 31.0 ms | 39.1 ms |
| SDNN | 44.4 ms | 60.3 ms |
| pNN50 | 9.3 % | 16.7 % |
| LF/HF | 1.15 | 1.15 |
| SD1 / SD2 | 22.0 / 58.6 ms | 27.7 / 80.8 ms |

兩日皆為規則竇性節律，指標落在健康成人短時記錄參考區間。8/20 資料在 t≈80–82 s 有一段動作雜訊，經 artifact 校正後不影響指標（詳見 commit 紀錄）。

> ⚠️ 本專案為課程教學用途。單導程、未校準振幅的自製裝置不具診斷效力，分析結果不構成醫療建議。
