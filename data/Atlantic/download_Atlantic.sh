#!/bin/bash

# Set the directory to save downloaded files
DOWNLOAD_DIR="sra"
FASTQ_DIR="fastq"
BGZIP_DIR="fastqgz"

# Create directories if they do not exist
mkdir -p $DOWNLOAD_DIR $FASTQ_DIR $BGZIP_DIR

# List of SRR files to download (example SRR IDs)
SRR_LIST="SRR_Acc_List_Atlantic.txt"

# Download SRR files from NCBI using SRA Toolkit and save them in DOWNLOAD_DIR
while read -r SRR; do
    echo "Downloading $SRR..."
    prefetch --max-size 80G "$SRR" -O "$DOWNLOAD_DIR"
done < "$SRR_LIST"

echo "Download completed."

# Convert downloaded SRR files to fastq and save them in FASTQ_DIR
for sra_file in $DOWNLOAD_DIR/*/*.sra; do
    echo "Converting $sra_file to fastq..."
    fastq-dump $sra_file --split-files -O $FASTQ_DIR/
done

echo "Fastq conversion completed."


# Convert fastq files to bgzip and save them in BGZIP_DIR
echo "Converting FASTQ files to BGZIP..."
for fastq_file in $FASTQ_DIR/*.fastq; do
  filename=$(basename "$fastq_file" .fastq)
  bgzip_file="$BGZIP_DIR/${filename}.fastq.gz"
  echo "Compressing $fastq_file to $bgzip_file"
  bgzip -c $fastq_file > $bgzip_file
done

echo "Conversion to bgzip completed."