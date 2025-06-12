#!/bin/bash

# fastqファイルをbgzipファイルを保存すディレクトリを指定

# fastq ファイルを保存するディレクトリを指定
FASTQ_DIR="fastq"

# ダウンロードしたいSRRのリストを指定します
SRR_LIST="SRR_Acc_list_test.txt"

# ディレクトリが存在しない場合は作成
mkdir -p $BGZIP_DIR

echo "Converting fastq files to bgzip..."
for fastq_file in $FASTQ_DIR/*.fastq; do
    filename=$(basename -- "$fastq_file")
    bgzip_file="$BGZIP_DIR/$filename.gz"
    echo "Converting $fastq_file to $bgzip_file"
    bgzip -c $fastq_file > $bgzip_file
done

echo "Conversion to bgzip completed."