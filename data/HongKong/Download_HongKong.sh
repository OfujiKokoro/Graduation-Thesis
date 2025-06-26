#!/bin/bash

# Set the directory to save downloaded files
DOWNLOAD_DIR="sra"
FASTQ_DIR="fastq"
BGZIP_DIR="fastqgz"
MERGE_FASTQ_DIR="merge_fastqgz"

# Create directories if they do not exist
mkdir -p $DOWNLOAD_DIR $FASTQ_DIR $BGZIP_DIR $MERGE_FASTQ_DIR

# List of SRR files to download
SRR_LIST="SRR_Acc_List_HongKong.txt"


# Download SRR files from NCBI using SRA Toolkit and save them in DOWNLOAD_DIR
while read -r SRR; do
    echo "Downloading $SRR..."
    prefetch --max-size 50G "$SRR" -O "$DOWNLOAD_DIR"
done < "$SRR_LIST"

echo "Download completed."

# Convert downloaded SRR files to fastq
for sra_file in $DOWNLOAD_DIR/*/*.sra; do
    echo "Converting $sra_file to fastq..."
    fastq-dump $sra_file --split-files -O $FASTQ_DIR/
done
echo "Fastq conversion completed."

# Convert fastq files to bgzip format
for fastq_file in $FASTQ_DIR/*.fastq; do
  filename=$(basename "$fastq_file" .fastq)
  bgzip_file="$BGZIP_DIR/${filename}.fastq.gz"
  echo "Compressing $fastq_file to $bgzip_file"
  bgzip -c $fastq_file > $bgzip_file
done

echo "Compression to bgzip completed."

# 1. Create list of fastq.gz files for _1 reads and concatenate
find "$BGZIP_DIR" -name '*_1.fastq.gz' | sort > "$MERGE_FASTQ_DIR/list_fastq_1_files.txt"
echo "Merging all _1.fastq.gz files into SRR446045_merge_1.fastq.gz..."
zcat $(cat "$MERGE_FASTQ_DIR/list_fastq_1_files.txt") | bgzip > "$MERGE_FASTQ_DIR/SRR446045_merge_1.fastq.gz"

# 2. Create list of fastq.gz files for _2 reads and concatenate
find "$BGZIP_DIR" -name '*_1.fastq.gz' | sort > "$MERGE_FASTQ_DIR/list_fastq_2s_files.txt"
echo "Merging all _1.fastq.gz files into SRR446045_merge_2.fastq.gz..."
zcat $(cat "$MERGE_FASTQ_DIR/list_fastq_1_files.txt") | bgzip > "$MERGE_FASTQ_DIR/SRR446045_merge_2.fastq.gz"

