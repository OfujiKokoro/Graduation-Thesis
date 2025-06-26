#!/bin/bash
set -eu

# Path to the reference genome file
reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

# Specify the directory where BAM files are stored
BAM_DIR=../02Mapping/bwa

# Specify the directory to save the variant calling results
VCFGZ_DIR="vcfgz"

mkdir -p "$VCFGZ_DIR"

CORE=5
min_MQ=20
min_BQ=20

# Perform variant calling for each BAM file
for bam_file in "$BAM_DIR"/*.bam; do
    # Extract sample name from BAM file name
    sample_name=$(basename "$bam_file" .bam)

    echo "Calling variants for $sample_name..."
    bcftools mpileup --threads $CORE --min-MQ $min_MQ --min-BQ $min_BQ \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    -f "$reference" "$bam_file" | \
    bcftools call --threads $CORE --format-fields GQ -m -Oz -o "$VCFGZ_DIR/$sample_name.vcf.gz"

    # Index the output VCF file
    bcftools index "$VCFGZ_DIR/$sample_name.vcf.gz"
done
