#!/bin/bash
set -euo pipefail

# Directory where FASTQ files are stored
TRIMMED_PAIRED_DIR=../01Qualiry_Control/trimmomatic/paired


# Directories to save results
PICARD_DIR="picard"
TMP_DIR="tmp"
FASTQ_DIR="fastq"
BGZIP_DIR="fastqgz"

# Create directories if they don't exist
mkdir -p "$PICARD_DIR" "$FASTQ_DIR" "$TMP_DIR" "$BGZIP_DIR"

# Get sample name from the first FASTQ file
fastq_file_1="$TRIMMED_PAIRED_DIR/SRR446045_merge_1.paired.fastq.gz"
sample_name=$(basename "$fastq_file_1" _1.paired.fastq.gz)
fastq_file_2="$TRIMMED_PAIRED_DIR/${sample_name}_2.paired.fastq.gz"

# Check if the paired file exists
if [[ ! -f "$fastq_file_2" ]]; then
    echo "Warning: Paired file for $fastq_file_1 not found. Exiting..."
    exit 1
fi

# Run Picard FastqToSam command
picard FastqToSam \
    F1="$fastq_file_1" \
    F2="$fastq_file_2" \
    O="$PICARD_DIR/${sample_name}_unaligned_reads.bam" \
    SM="$sample_name" \
    RG="$sample_name" \
    TMP_DIR="$TMP_DIR"

# Run bedtools bamtofastq command
bedtools bamtofastq -i "$PICARD_DIR/${sample_name}_unaligned_reads.bam" \
    -fq "$FASTQ_DIR/${sample_name}_1.fastq" \
    -fq2 "$FASTQ_DIR/${sample_name}_2.fastq"

# Convert fastq files to bgzip and save them in BGZIP_DIR
for fastq_file in $FASTQ_DIR/*.fastq; do
    filename=$(basename "$fastq_file" .fastq)
    bgzip_file="$BGZIP_DIR/$filename.fastq.gz"
    echo "Converting $fastq_file to $bgzip_file"
    bgzip -c $fastq_file > $bgzip_file
done
