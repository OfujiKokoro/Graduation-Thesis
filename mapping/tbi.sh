#!/bin/bash

# バリアントコール結果を保存するディレクトリを指定
VCFGZ_DIR="vcfgz"
VCFGZTBI_DIR="vcfgztbi"

# ディレクトリが存在しない場合は作成します
mkdir -p "$VCFGZ_DIR"
mkdir -p "$VCFGZTBI_DIR"

# BAM ファイルごとにバリアントをコール
for vcfgz_file in "$VCFGZ_DIR"/*.vcf.gz; do
    # サンプル名を取得
    sample_name=$(basename "$vcfgz_file" .vcf.gz)
    
    echo "Calling variants for $sample_name..."
    tabix -p vcf "$vcfgz_file"
    mv "$VCFGZ_DIR/$sample_name.vcf.gz.tbi" "$VCFGZTBI_DIR"
done

echo  "Variant calling completed."
