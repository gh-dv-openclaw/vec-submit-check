#!/usr/bin/env bash
# 拉取每块板的基因 panel 和机读索引。这些是公开的，不需要登录。
#
# panels_index.json 里除了基因表还有每块板的细胞数区间、必需的 obsm 键，
# 以及每个指标的 floor/ceiling 锚点。
set -euo pipefail
cd "$(dirname "$0")"

SITE=https://virtualembryo.ai
mkdir -p panels

curl -fsS --max-time 60 "$SITE/challenge/panels/index.json" -o panels/panels_index.json

for f in $(python3 -c '
import json
print(" ".join(v["genes_file"] for v in json.load(open("panels/panels_index.json")).values()))
'); do
  curl -fsS --max-time 120 "$SITE/challenge/panels/$f" -o "panels/$f"
  echo "panels/$f: $(wc -l < "panels/$f") genes"
done
