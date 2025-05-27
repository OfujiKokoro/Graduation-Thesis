#!/bin/bash
set -eu

# バリアントコール結果を保存するディレクトリを指定b
BCF_DIR="bcf"
VCF_SNP_DIR="vcf_snp"
VCFGZ_SNP_DIR="vcfgz_snp"

# ディレクトリが存在しない場合は作成します
mkdir -p "$VCF_SNP_DIR"
mkdir -p "$VCFGZ_SNP_DIR"

# BAM ファイルごとにバリアントをコール
for bcf_file in "$BCF_DIR"/*.bcf; do
    # サンプル名を取得
    sample_name=$(basename "$bcf_file" .bcf)
    
    echo "Calling variants for $sample_name..."
    bcftools call --threads 6 -vm -Ov -o "$VCF_SNP_DIR/$sample_name.vcf" "$BCF_DIR/$sample_name.bcf"
    bcftools call --threads 6 -vm -Oz -o "$VCFGZ_SNP_DIR/$sample_name.vcf.gz" "$BCF_DIR/$sample_name.bcf"
    bcftools index "$VCFGZ_SNP_DIR/$sample_name.vcf.gz"
done

echo  "Variant calling completed."
