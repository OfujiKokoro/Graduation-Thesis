#!/bin/bash
set -eu

# バリアントコール結果を保存するディレクトリを指定
VCF_DIR="vcf"
BED_VCF_DIR="bed_vcf"  # ディレクトリ名を修正

# ディレクトリが存在しない場合は作成します
mkdir -p "$VCF_DIR"
mkdir -p "$BED_VCF_DIR"

# BAM ファイルごとにバリアントをコール
for vcf_file in "$VCF_DIR"/*.vcf; do
    # サンプル名を取得
    sample_name=$(basename "$vcf_file" .vcf)
    
    echo "Filtering vcf file for $sample_name..."
    cat <(grep '^#' "$vcf_file") <(bedtools intersect -a "$vcf_file" -b filtered_100_samples_mini_position.bed) > "$BED_VCF_DIR/$sample_name.vcf"
    echo "Filtering completed for $sample_name" 
done

echo "Variant calling completed."
