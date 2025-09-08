#!/bin/bash
set -eu

# Specify the directory to store BAM files after MarkDuplicates
MARK_BAM_DIR="markdup_bam"
METRICS_DIR="metrics_txt"
BAM_DP_DIR="bam_DP"
mkdir -p "$BAM_DP_DIR" "$MARK_BAM_DIR" "$METRICS_DIR"

reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

CORE=5

process_dp_vcf() {
    bam_file="$1"
    sample_name=$(basename "$bam_file" .markdup.bam)
    echo "Processing BAM file for sample: $sample_name"
    BAM_DP_DIR="bam_DP"

    reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

    picard -Xms40G -Xmx50G CollectWgsMetrics \
    I="$bam_file" \
    O="$BAM_DP_DIR/${sample_name}_depth_metrics.txt" \
    R="$reference"
}

export -f process_dp_vcf

# Use find to locate BAM files and parallel to execute them in parallel (maximum of 2 jobs simultaneously)
find "$MARK_BAM_DIR" -name "*.bam" | parallel -j 5 process_dp_vcf
echo "All BAM files have been processed."
