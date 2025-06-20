#!/bin/bash

# Set directories for downloaded and processed files
DOWNLOAD_DIR="sra"
FASTQ_DIR="fastq"
BGZIP_DIR="fastqgz"

# Create directories if they do not exist
mkdir -p $DOWNLOAD_DIR $FASTQ_DIR $BGZIP_DIR

# File containing the list of SRR accessions to download
SRR_LIST="SRR_Acc_List_Pacific.txt"

# Download SRR files from NCBI using SRA Toolkit and save them in DOWNLOAD_DIR
while read -r SRR; do
  echo "Downloading $SRR..."
  prefetch --max-size 20G "$SRR" -O "$DOWNLOAD_DIR"
done < "$SRR_LIST"

# Convert downloaded .sra files to fastq and save in the fastq directory
for file in $DOWNLOAD_DIR/*/*.sra; do
  echo "Converting $file to FASTQ..."
  fastq-dump $file --split-files -O $FASTQ_DIR/
  echo "FASTQ conversion completed for $file."
done

echo "Download and FASTQ conversion completed."

# File containing list of SRR replacements (old\tnew)
REPLACE_LIST="Replace_List.txt"

# Define function to replace SRR IDs in filenames and content
Replace_SRR() {
  line="$1"
  old_string=$(echo "$line" | cut -f1)
  new_string=$(echo "$line" | cut -f2)

  input_file="$FASTQ_DIR/${old_string}_2.fastq"
  output_file="$FASTQ_DIR/${new_string}_2.fastq"

  # Replace SRR ID in file contents in-place
  sed -i "s/$old_string/$new_string/g" "$input_file" 
  mv "$input_file" "$output_file"
  
  echo "Replacement completed: $old_string -> $new_string"
}

# Read the replacement list and perform replacements
while IFS= read -r line; do
  Replace_SRR "$line"
done < "$REPLACE_LIST"

# Convert all fastq files to bgzip-compressed files and save in the bgzip directory
echo "Converting FASTQ files to BGZIP..."
for fastq_file in $FASTQ_DIR/*.fastq; do
  filename=$(basename -- "$fastq_file")
  filename_without_extension="${filename%.*}"
  bgzip_file="$BGZIP_DIR/$filename_without_extension.fastq.gz"
  echo "Compressing $fastq_file to $bgzip_file"
  bgzip -c $fastq_file > $bgzip_file
done

echo "BGZIP compression completed."
