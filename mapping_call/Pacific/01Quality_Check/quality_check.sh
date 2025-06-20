#!/bin/bash
set -eu  # Enable error handling: exit on error or undefined variable

# Define the directory containing FASTQ files (relative path)
PACIFIC_DATA_DIR=../../../data/Pacific/fastqgz

# Define output directories for FastQC and MultiQC results
QC_DIR="quality_result"
MULTI_QC_DIR="multiqc"

# Create output directories if they do not exist
mkdir -p "$QC_DIR" "$MULTI_QC_DIR"

# Run FastQC on each .fastq.gz file in the input directory
for fastqgz_file in "$PACIFIC_DATA_DIR"/*.fastq.gz; do
    fastqc "$fastqgz_file" -o "$QC_DIR"
done

# Run MultiQC to aggregate FastQC reports
multiqc "$QC_DIR" -o "$MULTI_QC_DIR"
