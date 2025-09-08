#!/bin/bash
#GATK3
set -euo pipefail

BAM_DIR=../04MarkDuplicate/markdup_bam
GATK_DIR="."
INTERVAL_DIR="$GATK_DIR/interval"
BAM_OUT_DIR="$GATK_DIR/bam_out"
GVCF="$GATK_DIR/GenotypeGVCFs"
VCF="$GATK_DIR/VCF"
TVCF="$GATK_DIR/SelectVariants"
MYLOG="$GATK_DIR/log"

mkdir -p "$INTERVAL_DIR" "$BAM_OUT_DIR" "$GVCF" "$VCF" "$TVCF" "$MYLOG"


# Path to the reference genome file
reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

#conda activate gatk4

#gatk CreateSequenceDictionary \
#-R $reference -O ../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.dict

#conda activate gatk
process_vcf() {
    bam_file="$1"
    sample_name=$(basename "$bam_file" .markdup.bam)
    echo "Processing BAM file for sample: $sample_name"

    GATK_DIR="."
    INTERVAL_DIR="$GATK_DIR/interval"
    BAM_OUT_DIR="$GATK_DIR/bam_out"
    reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

    gatk3 -Xms30G -Xmx50G -T RealignerTargetCreator \
    -R $reference \
    -I $bam_file \
    -o $INTERVAL_DIR/${sample_name}.markdup.intervals \
    -drf BadMate
    
    gatk3 -Xms30G -Xmx50G -T IndelRealigner \
    -R $reference \
    -I $bam_file \
    -targetIntervals $INTERVAL_DIR/${sample_name}.markdup.intervals \
    --consensusDeterminationModel USE_READS  \
    -o $BAM_OUT_DIR/${sample_name}_markdup_idnel_realign.bam
    
    echo "Markduplicate $sample_name completed."
}

export -f process_vcf

# Use find to locate BAM files and parallel to execute them in parallel (maximum of 2 jobs simultaneously)
echo "Starting processing of BAM files in directory: $BAM_DIR..."

find "$BAM_DIR" -name "*.bam" | parallel -j 5 process_vcf

