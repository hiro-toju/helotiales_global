# 0_Data_Property

このステップでは、`4_Specificity_Lat` の前処理で除外された unique occurrence と植物ラベルを反映し、解析対象データの基本的な性質を確認するための表と図を作成します。

## 入力

通常は `4_Specificity_Lat` の出力runを指定します。

- `4_Specificity_Lat/outputs/<run_name>/data/helotiales_root_occurrences_with_hosts.tsv`
- `4_Specificity_Lat/outputs/<run_name>/data/plant_taxonomy.tsv`
- `4_Specificity_Lat/outputs/<run_name>/helotiales_root_occurrences_with_hosts/helotiales_root_occurrences_with_hosts_excluded_unique_occurrences.tsv`
- `4_Specificity_Lat/outputs/<run_name>/helotiales_root_occurrences_with_hosts/helotiales_root_occurrences_with_hosts_excluded_plant_columns.tsv`

## 基本コマンド

SH版の例です。

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 data_property_summary.py \
  --specificity-run-name sh150_family30_10000_lat20_minblock100 \
  --plot-top-n 15 \
  --top-sh-bar-n 50 \
  --figure-font-scale 1.2
```

真菌属版の例です。

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 data_property_summary.py \
  --specificity-run-name genus150_family30_10000_lat20_minblock100 \
  --plot-top-n 15 \
  --top-sh-bar-n 50 \
  --figure-font-scale 1.2
```

`data_property_summary.py` はrun名で入力フォルダを切り替えます。そのため、SH版と真菌属版で使うPythonコードは同じです。

## 主なオプション

| オプション | 内容 |
|---|---|
| `--specificity-run-name` | 参照する `4_Specificity_Lat` のrun名 |
| `--input` | occurrence TSVを明示指定する場合に使用 |
| `--specificity-result-dir` | 除外レポートがあるフォルダを明示指定する場合に使用 |
| `--plant-taxonomy` | 植物分類表を明示指定する場合に使用 |
| `--output-dir` | 出力先。省略時は `outputs/<run_name>` |
| `--plot-top-n` | 円グラフ等で個別表示する上位カテゴリ数 |
| `--top-sh-bar-n` | SH出現頻度の棒グラフに表示する上位SH数。デフォルト50 |
| `--figure-font-scale` | 図中文字サイズの倍率 |

## 出力

主な出力は `outputs/<run_name>/` に作成されます。

- 除外後のcleaned occurrence table
- fungal SH, fungal genus, fungal family, plant genus, plant family, plant order, country, continent の集計表
- 各カテゴリの unique occurrence 内訳を示すPDF図
- 出現頻度上位SHの棒グラフPDF
- 実行コマンド、計算時間、Python・パッケージのバージョン情報

## 注意

このステップは、下流解析の入力を作り替えるための本体ステップではなく、データの性質を確認するための補助的な集計ステップです。
