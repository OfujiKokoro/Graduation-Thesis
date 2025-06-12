#!/bin/bash

# 入力FASTQファイル名
input_file="/Users/kokoro/lab/seeturtle/Test4/DRR217740_2.fastq"

# 置換後のFASTQファイル名
output_file="DRR2177439_2.fastq"

# 置換する文字列
old_string="DRR217740"
new_string="DRR217739"

# sedコマンドを使って置換を実行し、結果をoutput.fastqに書き込む
sed "s/$old_string/$new_string/g" "$input_file" > "$output_file"

echo "置換が完了しました。"
