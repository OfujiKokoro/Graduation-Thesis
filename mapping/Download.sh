#!/bin/bash

# ダウンロードしたいSRRのリストを指定
SRR_LIST="SRR_Acc_List.txt"

# ダウンロードしたファイルを保存するディレクトリを指定
DOWNLOAD_DIR="downloads"

# bgzipファイルを保存するディレクトリを指定
BGZIP_DIR="bgzip"

# fastq ファイルを保存するディレクトリを指定
FASTQ_DIR="fastq"

# ディレクトリが存在しない場合は作成
mkdir -p $DOWNLOAD_DIR
mkdir -p $FASTQ_DIR
mkdir -p $BGZIP_DIR

# SRR リストから一括でダウンロード
echo "Downloading SRR files..."
cat $SRR_LIST | xargs prefetch -O $DOWNLOAD_DIR

# ダウンロードしたSRRファイルから fastq ファイルを作成、fastq ディレクトリに保存
for file in $DOWNLOAD_DIR/*/*.sra; do
    echo "Converting $file to fastq..."
    fastq-dump $file --split-files -O $FASTQ_DIR/
	echo "Converting fastq files to bgzip..."

done

echo "Download and fastq conversion completed."

# 作成したfastq ファイルを　bgzipファイルに変換して、bzgipディレクトリに保存
echo "Converting fastq files to bgzip..."
for fastq_file in $FASTQ_DIR/*.fastq; do
    filename=$(basename -- "$fastq_file")
    filename_without_extension="${filename%.*}"
    bgzip_file="$BGZIP_DIR/$filename_without_extension.fastq.gz"
    echo "Converting $fastq_file to $bgzip_file"
    bgzip -c $fastq_file > $bgzip_file
done

echo "Conversion to bgzip completed."
