#!/bin/bash

# バリアントコール結果を保存するディレクトリを指定
VCF_DIR="vcf_filter_DP>=3"
VCFDP_DIR="vcf_query_DP3"

# ディレクトリが存在しない場合は作成します
mkdir -p "$VCF_DIR"
mkdir -p "$VCFDP_DIR"  # ディレクトリ名が誤っていたため修正

for vcf_file in "$VCF_DIR"/*.vcf; do  # vcf_file の変数名が正しくないため修正
    # サンプル名を取得
    sample_name=$(basename "$vcf_file" .vcf)
    echo "Query vcffile for $sample_name..."
    # 出力ファイル名を定義します
    output_file="${VCFDP_DIR}/${sample_name}_DP3.txt"  # 出力ファイル名の変数が正しくないため修正
    bcftools query -f '%DP\n' "$vcf_file" > "$output_file"
    echo "DP の抽出が完了しました: $output_file"
done

echo  "Variant calling completed."
