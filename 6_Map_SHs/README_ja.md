# 6_Map_SHs

`1_GlobalFungi` で作成した root unique occurrence の表から、出現数の多いSHを選び、SHごとの地理分布図をPDFで出力する独立ツールです。

旧 `2_Specificity` に含まれていた `*_selected_SH_spatial_maps/` の機能だけを切り出したものです。randomization、d′、2DP、FASTA抽出は行いません。

## 基本コマンド

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 map_abundant_shs.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --sh-selection rank \
  --sh-rank-max 50 \
  --figure-font-scale 1.5
```

この例では `outputs/sh50_maps/selected_SH_spatial_maps/` に、上位50 SHのPDF地図を出力します。

## 主なオプション

- `--sh-selection rank`: unique occurrence数の順位でSHを選択します。
- `--sh-rank-max 50`: 上位何SHを出力するかを指定します。
- `--sh-selection min-count --sh-min-count 100`: unique occurrence数が指定値以上のSHをすべて出力します。
- `--sh-selection list --selected-sh-file selected_shs.tsv`: 明示的に指定したSHだけを出力します。
- `--keep-uncertain-hosts`: 曖昧な植物ラベルを含むoccurrenceも地図に含めます。指定しない場合は、旧 `2_Specificity` の地図に近い挙動として、不確かなhost labelを除外します。
- `--figure-font-scale`: 図中の文字サイズ倍率です。
- `--point-size`: 地図上の点の大きさです。

## 出力

- `selected_SH_spatial_maps/*.pdf`: SHごとの地理分布図
- `sh_unique_occurrence_ranking.tsv`: 全SHのunique occurrence数ランキング
- `selected_shs.tsv`: 地図を出力したSHの一覧
- `excluded_occurrences.tsv`: 地図から除外したoccurrence
- `run_summary.tsv`: 入力行数、除外行数、出力SH数、計算時間
- `command.txt`: 実行コマンド
- `package_versions.tsv`: Pythonと使用パッケージのバージョン
- `analysis.log`: 簡易ログ

## カウント単位

SHの順位付けは、`sh_id + latitude + longitude + host_candidate` で定義される unique occurrence 数に基づきます。これは、同じ緯度・経度・host candidateでまとめられたsampling unit相当の出現数です。

