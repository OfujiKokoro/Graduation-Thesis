#!/bin/bash
set -euo pipefail


accession=GCF_015237465.2
species=Chelonia_mydas



datasets download genome accession "$accession" \
  --include genome,gff3,cds,rna,protein,seq-report \
   --filename "${accession}.zip"

unzip -o "${accession}.zip"

mv ncbi_dataset/data ./ 
rm -r ncbi_dataset  

TARGET_DIR=data/${accession}
  
for file in "$TARGET_DIR"/*; do
  base=$(basename "$file")
  newname="${species}_${base}"
  mv "$file" "$TARGET_DIR/$newname"
done

