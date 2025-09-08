#!/bin/bash
# GATK4
set -euo pipefail

GATK_DIR="."
INTERVAL_DIR="$GATK_DIR/interval"
BAM_OUT_DIR="$GATK_DIR/bam_out"
GVCF="$GATK_DIR/GenotypeGVCFs"
VCF="$GATK_DIR/VCF"
TVCF="$GATK_DIR/SelectVariants"
MYLOG="$GATK_DIR/log"

mkdir -p "$INTERVAL_DIR" "$BAM_OUT_DIR" "$GVCF" "$VCF" "$TVCF" "$MYLOG"

# conda activate gatk4

reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

for bam_file in "$BAM_OUT_DIR"/*.bam; do
    sample_name=$(basename "$bam_file" _markdup_idnel_realign.bam)
    echo "Processing BAM file for sample: $sample_name"

    touch "$MYLOG/${sample_name}_gatk.log"

    gatk HaplotypeCaller \
        --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        -R "$reference" \
        -ERC BP_RESOLUTION \
        -mbq 20 \
        -I "$bam_file" \
        -O "$GVCF/${sample_name}.g.vcf.gz" \
        --output-mode EMIT_ALL_ACTIVE_SITES \
        >> "$MYLOG/${sample_name}_gatk.log" 2>&1

    gatk GenotypeGVCFs \
        --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        -R "$reference" \
        --include-non-variant-sites \
        --standard-min-confidence-threshold-for-calling 0 \
        -V "$GVCF/${sample_name}.g.vcf.gz" \
        -O "$VCF/${sample_name}.vcf.gz" \
        >> "$MYLOG/${sample_name}_gatk.log" 2>&1

    gatk SelectVariants \
        --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
        -R "$reference" \
        --remove-unused-alternates \
        -V "$VCF/${sample_name}.vcf.gz" \
        -O "$TVCF/${sample_name}.TrimAlt.vcf.gz" \
        >> "$MYLOG/${sample_name}_gatk.log" 2>&1

    echo "$sample_name completed."
done
