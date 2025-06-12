#!/bin/bash#!/bin/bash

# 出力ファイル
output_file="vcf_stats_summary.txt"

# VCFファイルのパス
vcfgz_dir="/home/ofuji/Lab/seaturtle/vcfgz/*.vcf.gz"

# 各VCFファイルについて統計情報を取得して出力ファイルに書き込む
for vcfgz_file in $vcfgz_dir; do
    file_name=$(basename "$vcfgz_file")
    stats=$(rtg vcfstats "$vcfgz_file")
    echo -e "$file_name\t$stats" >> "$output_file"
done

echo "Summary statistics written to $output_file"
