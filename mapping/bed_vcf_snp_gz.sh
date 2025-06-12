#!/bin/bash
set -eu

# fastq ファイルを保存するディレクトリを指定
BED_VCF_SNP_DIR="bed_vcf_snp"
BED_VCF_SNP_GZ_DIR="bed_vcfgz_snp"

# ディレクトリが存在しない場合は作成
mkdir -p "$BED_VCF_SNP_DIR"
mkdir -p "$BED_VCF_SNP_GZ_DIR"

echo "Converting to bgzip..."
for vcf_file in "$BED_VCF_SNP_DIR"/*.vcf; do
    filename=$(basename "$vcf_file" .vcf)
    echo "Converting $vcf_file to $vcf_file.gz"
    bgzip "$vcf_file"
done

echo "Conversion to bgzip completed."
