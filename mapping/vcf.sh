#!/bin/bash

# レファレンスファイルのパスを指定
REFERENCE="/home/ofuji/Lab/seaturtle/reference/ncbi_dataset/data/GCF_015237465.2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna"

# BAM ファイルが保存されているディレクトリを指定
BAM_DIR="bam"

# バリアントコール結果を保存するディレクトリを指定b
BCF_DIR="bcf"
VCF_DIR="vcf"
VCFGZ_DIR="vcfgz"

# ディレクトリが存在しない場合は作成します
mkdir -p "$BCF_DIR"
mkdir -p "$VCF_DIR"
mkdir -p "$VCFGZ_DIR"

# BAM ファイルごとにバリアントをコール
for bam_file in "$BAM_DIR"/*.bam; do
    # サンプル名を取得
    sample_name=$(basename "$bam_file" .bam)
    
    echo "Calling variants for $sample_name..."
    bcftools mpileup --threads 6 -Ob -o "$BCF_DIR/$sample_name.bcf" -f "$REFERENCE" "$bam_file" 
    bcftools call --threads 6 -m -Ov -o "$VCF_DIR/$sample_name.vcf" "$BCF_DIR/$sample_name.bcf"
    bcftools call --threads 6 -m -Oz -o "$VCFGZ_DIR/$sample_name.vcf.gz" "$BCF_DIR/$sample_name.bcf"
    bcftools index "$VCFGZ_DIR/$sample_name.vcf.gz"
done

echo  "Variant calling completed."
