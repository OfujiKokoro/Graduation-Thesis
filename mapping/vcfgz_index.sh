#!/bin/bash
set -eu

VCFGZ_DIR="bed_vcfgz_snp"

for vcf_file in "$VCFGZ_DIR"/*.vcf.gz; do
    filename=$(basename "$vcf_file" .vcf.gz)
    echo "Calling variants for $filename..."
    bcftools index "$vcf_file"
done

echo  "Variant calling completed."
