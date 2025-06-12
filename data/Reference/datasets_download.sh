#!/bin/bash
set -euo pipefail

# 並列実行数（コア数）
NUM_CORES=15

# アクセッションリストのファイルパス
SPECIES_LIST=/home/ofuji/Lab/Molecular_evolution/Datasets/Species_list.txt
ACESSION_LIST=/home/ofuji/Lab/Molecular_evolution/Datasets/Acession_Species_GCA_GCF_list.txt

while read -r species; do
  echo "Creating directory for $species"
  mkdir -p "$species"
done < "$SPECIES_LIST"


# ダウンロード処理の関数を定義
process_datasets() {
  line="$1"
  accession=$(echo "$line" | cut -f1)
  species=$(echo "$line" | cut -f2)

  cd ${species}
  echo "Downloading $accession ($species)..."
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
  
}


export -f process_datasets

# GNU parallel で並列実行
cat "$ACESSION_LIST" | parallel -j "$NUM_CORES" process_datasets {}
