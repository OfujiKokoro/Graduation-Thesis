#!/bin/bash

# バリアントコール結果を保存するディレクトリを指定
VCF_DIR="vcf"
VCFFILTER_DIR="vcf_filter"

# ディレクトリが存在しない場合は作成します
mkdir -p "$VCF_DIR"
mkdir -p "$VCFFILTER_DIR"  # ディレクトリ名が誤っていたため修正

# BAM ファイルごとにバリアントをコール
for vcf_file in "$VCF_DIR"/*.vcf; do  # vcf_file の変数名が正しくないため修正
    # サンプル名を取得
    sample_name=$(basename "$vcf_file" .vcf)
    
    echo "Filtering vcffile for $sample_name..."
    bcftools filter -i 'DP>=3' -o "$VCFFILTER_DIR/$sample_name.vcf" "$vcf_file"  # 出力ファイル名に修正
    
    echo "Filtering completed for $sample_name"  # フィルタリングが完了したことを表示
done

echo  "Variant calling completed."
