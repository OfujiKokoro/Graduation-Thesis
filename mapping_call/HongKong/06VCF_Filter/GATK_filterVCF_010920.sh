#!/bin/bash
set -eu

# Directories for saving variant call results
RENAME_DIR="01rename"
INDEL_DIR="02Indel"
DP_GQ_DIR="03DP_GQ"

# Create directories if they don't exist
mkdir -p "$RENAME_DIR" "$INDEL_DIR" "$DP_GQ_DIR"

VCF_DIR="../05GATK_Call/SelectVariants"
CHR_LIST="rename_chr_noX.txt"
METRICS_DIR="../04MarkDuplicate/bam_DP"
python_script="filterVCF_010920.py"

for vcf_file in "$VCF_DIR"/*.vcf.gz; do
  sample_name=$(basename "$vcf_file" _merge.TrimAlt.vcf.gz)
  metrics_file="$METRICS_DIR/${sample_name}_merge_depth_metrics.txt"

  # Step 1: Rename chromosome names according to the provided list
  bcftools annotate --rename-chrs "$CHR_LIST" "$vcf_file"  \
  -Oz -o "${RENAME_DIR}/${sample_name}_renamed.vcf.gz"
  
  # Step 2: Remove indels, keep only SNPs
  vcftools --gzvcf "${RENAME_DIR}/${sample_name}_renamed.vcf.gz" --remove-indels \
  --recode --recode-INFO-all \
  --stdout | bgzip > "${INDEL_DIR}/${sample_name}_rm_indel.vcf.gz"

  # Extract MEAN_COVERAGE from the metrics file
  mean_cov=$(awk '/^## METRICS CLASS/{getline; getline; print $2}' "$metrics_file")

  # Define DP thresholds (1/3 × mean coverage ~ 2 × mean coverage)
  MIN_DP=$(awk -v mc="$mean_cov" 'BEGIN{printf "%.0f", mc/3}')
  MAX_DP=$(awk -v mc="$mean_cov" 'BEGIN{printf "%.0f", mc*2}')

  echo "$sample_name : MEAN_COVERAGE = $mean_cov (DP range = $MIN_DP - $MAX_DP)"

  # Minimum genotype quality (GQ) threshold
  MIN_GQ=20

  # Step 3: Apply DP and GQ filtering using the custom Python script
  python3 "$python_script" "${INDEL_DIR}/${sample_name}_rm_indel.vcf.gz" \
  "$MIN_DP" "$MAX_DP" | bgzip > "${DP_GQ_DIR}/${sample_name}_rm_indel_DP${MIN_DP}_${MAX_DP}_GQ${MIN_GQ}.vcf.gz"

  # Index the filtered VCF
  tabix -p vcf "${DP_GQ_DIR}/${sample_name}_rm_indel_DP${MIN_DP}_${MAX_DP}_GQ${MIN_GQ}.vcf.gz"
done

echo "All samples have been processed."
