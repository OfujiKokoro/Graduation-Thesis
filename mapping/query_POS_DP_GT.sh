#!/bin/bash

# バリアントコール結果を保存するディレクトリを指定
VCF_DIR="vcf"
VCFDP_DIR="vcf_query_POS_DP_GT"

# ディレクトリが存在しない場合は作成します
mkdir -p "$VCF_DIR"
mkdir -p "$VCFDP_DIR"

for vcf_file in "$VCF_DIR"/*.vcf; do
    # サンプル名を取得
    sample_name=$(basename "$vcf_file" .vcf)
    echo "Query vcffile for $sample_name..."
    # 出力ファイル名を定義します
    output_file="${VCFDP_DIR}/${sample_name}_POS_DP_GT.txt"
    # bcftools query コマンドを修正し、ファイルパスを指定して実行します
    bcftools query -f '%CHROM %POS %DP[\t%GT]\n' "$vcf_file" > "$output_file"

    echo "DP の抽出が完了しました: $output_file"
done

echo "Variant calling completed."
