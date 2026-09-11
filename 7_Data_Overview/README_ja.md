# 7_Data_Overview

Helotiales occurrence データの概要図を作成する補助ワークフローです。

以下の4つのデータ定義について、サンプリング地点の世界地図と continent 構成比のPDF図を出力します。

1. GlobalFungi の全 occurrence
2. Helotiales が含まれる全 occurrence
3. 上記のうち、`sample_type == root` の全 occurrence
4. 信頼できる宿主植物情報がある root unique occurrence

1番目の GlobalFungi 全 occurrence は、`GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz` 内の sample × SH abundance matrix において、値が正のセルを1 occurrenceとして数えます。地図では、少なくとも1つのSHが検出された sample の緯度・経度をプロットします。

4番目のデータは、既定では `0_Data_Property` で198 unique occurrencesを除外した cleaned table を使用します。該当ファイルがない場合は、`1_GlobalFungi/results_helotiales/helotiales_root_occurrences_with_hosts.tsv` を使用します。

## 基本コマンド

```bash
cd "[CURRENT DIRECTORY]"
python3 data_overview.py \
  --globalfungi-data ../1_GlobalFungi/data \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --figure-font-scale 1.5
```

## 主な出力

- `outputs/data_overview_summary.tsv`  
  各データ定義の occurrence 数、座標つき occurrence 数、サンプリング地点数、continent 数。

- `outputs/continent_breakdown.tsv`  
  各データ定義における continent ごとの count と proportion。

- `outputs/figures/*_world_map.pdf`  
  各データ定義のサンプリング地点マップ。

- `outputs/figures/*_continent_proportions.pdf`  
  各データ定義の continent 構成比。

- `outputs/run_configuration.tsv`, `runtime_versions.tsv`, `run_timing.tsv`, `workflow.log`  
  再現性確認用の設定・バージョン・計算時間・ログ。
