#!/bin/bash
set -eu

# Specify the directory where BAM files are stored
BAM_DIR=../03Mapping/bwa
# Specify the directories to store outputs
MARK_BAM_DIR="markdup_bam"
METRICS_DIR="metrics_txt"
BAM_DP_DIR="bam_DP"
mkdir -p "$BAM_DP_DIR" "$MARK_BAM_DIR" "$METRICS_DIR"

reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

CORE=5

echo "Starting processing of BAM files in directory: $BAM_DIR..."

for bam_file in "$BAM_DIR"/*.bam; do
    [ ! -f "$bam_file" ] && echo "No BAM files found in $BAM_DIR." && break

    sample_name=$(basename "$bam_file" .bam)
    echo "Processing BAM file for sample: $sample_name"

    export 
    TMPDIR=/home/ofuji/github_test/mapping_call/Atlantic/03MarkDuplicate/tmp 
    mkdir -p "$TMPDIR"


    # Step 1: MarkDuplicates
    picard -Xms100G -Xmx150G -XX:ParallelGCThreads=$CORE MarkDuplicates \
        -Djava.io.tmpdir=$TMPDIR \
        I="$bam_file" \
        O="$MARK_BAM_DIR/${sample_name}.markdup.bam" \
        M="$METRICS_DIR/${sample_name}.markdup.metrics.txt" \
        ASSUME_SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=TRUE \
        MAX_RECORDS_IN_RAM=50000000

    samtools index -@ "$CORE" \
        "$MARK_BAM_DIR/${sample_name}.markdup.bam" \
        "$MARK_BAM_DIR/${sample_name}.markdup.bam.bai"

    # Step 2: CollectWgsMetrics
    picard -Xms100G -Xmx150G CollectWgsMetrics \
        I="$MARK_BAM_DIR/${sample_name}.markdup.bam" \
        O="$BAM_DP_DIR/${sample_name}_depth_metrics.txt" \
        R="$reference"

    rm -rf "$TMPDIR"/*
done