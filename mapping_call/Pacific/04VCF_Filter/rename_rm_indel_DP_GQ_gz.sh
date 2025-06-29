#!/bin/bash
set -eu

# Directories for saving variant call results
RENAME_DIR="01rename"
INDEL_DIR="02Indel"
DP_DIR="03DP"
GQ_DIR="04GQ"

# Create directories if they don't exist
mkdir -p "$RENAME_DIR" "$INDEL_DIR" "$DP_DIR" "$GQ_DIR"

# Number of cores
CORE=16

# DP filter thresholds
MIN_DP=10
MAX_DP=200
MIN_GQ=30


# Original file path
VCF_DIR="../03Bcftools_Call/vcfgz"
CHR_LIST="rename_chr_noX.txt"

# Process each VCF file
for vcf_file in "$VCF_DIR"/*.vcf.gz; do
  sample_name=$(basename "$vcf_file" .vcf.gz)
  
  echo "Processing sample: $sample_name"

  # Step 1: Rename chromosomes
  bcftools annotate --rename-chrs "$CHR_LIST" "$vcf_file" --threads "$CORE" \
  -Oz -o "${RENAME_DIR}/${sample_name}_renamed.vcf.gz"
  
  # Step 2: Remove indels
  vcftools --gzvcf "${RENAME_DIR}/${sample_name}_renamed.vcf.gz" --remove-indels --recode --recode-INFO-all \
  --stdout | bgzip > "${INDEL_DIR}/${sample_name}_rm_indel.vcf.gz"
  
  # Step 3: DP filter
  bcftools filter -i "FORMAT/DP>=$MIN_DP & FORMAT/DP<=$MAX_DP" "${INDEL_DIR}/${sample_name}_rm_indel.vcf.gz" \
  -Oz -o "${DP_DIR}/${sample_name}_rm_indel_DP${MIN_DP}_${MAX_DP}.vcf.gz" 

  # Step 4: GQ filter
  bcftools filter -e "FORMAT/GQ<$MIN_GQ" "${DP_DIR}/${sample_name}_rm_indel_DP${MIN_DP}_${MAX_DP}.vcf.gz" \
  -Oz -o "${GQ_DIR}/${sample_name}_rm_indel_DP${MIN_DP}_${MAX_DP}_GQ${MIN_GQ}.vcf.gz" 
done

echo "All samples processed successfully."



