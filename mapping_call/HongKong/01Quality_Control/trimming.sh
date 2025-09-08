
#!/bin/bash
set -eu  

# Define directory paths
HONGKONG_DATA_DIR=../../../data/HongKong/merge_fastqgz
TRIMMED_DIR="trimmomatic"
TRIMMED_PAIRED_DIR="$TRIMMED_DIR/paired"
TRIMMED_UNPAIRED_DIR="$TRIMMED_DIR/unpaired"
QC_DIR="quality_result"


CORE=10  

# 必要なディレクトリが存在しない場合は作成
mkdir -p "$TRIMMED_DIR" "$TRIMMED_PAIRED_DIR" "$TRIMMED_UNPAIRED_DIR" "$QC_DIR"

# Trimming HongKong FASTQ files using Trimmomatic
for fastq_file_1 in "$HONGKONG_DATA_DIR"/*_1.fastq.gz; do
    sample_name=$(basename "$fastq_file_1" _1.fastq.gz)
    fastq_file_2="$HONGKONG_DATA_DIR/${sample_name}_2.fastq.gz"
    
    if [[ ! -f "$fastq_file_2" ]]; then
        echo "Warning: Paired file for $fastq_file_1 not found. Skipping..."
        continue
    fi

    echo "Trimming $sample_name start..."

    trimmed_fastq_file_1="$TRIMMED_PAIRED_DIR/${sample_name}_1.paired.fastq.gz"
    trimmed_fastq_file_2="$TRIMMED_PAIRED_DIR/${sample_name}_2.paired.fastq.gz"
    unpaired_fastq_file_1="$TRIMMED_UNPAIRED_DIR/${sample_name}_1.unpaired.fastq.gz"
    unpaired_fastq_file_2="$TRIMMED_UNPAIRED_DIR/${sample_name}_2.unpaired.fastq.gz"

    trimmomatic PE -phred33 -threads $CORE "$fastq_file_1" "$fastq_file_2" \
        "$trimmed_fastq_file_1" "$unpaired_fastq_file_1" \
        "$trimmed_fastq_file_2" "$unpaired_fastq_file_2" \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
        LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:30

    echo "Trimming $sample_name completed."
done

echo "Trimming all samples completed."

# Run FastQC on each trimmed paired .fastq.gz file
for fastqgz_file in "$TRIMMED_PAIRED_DIR"/*.fastq.gz; do
    echo "Running FastQC on $fastqgz_file..."
    fastqc "$fastqgz_file" -o "$QC_DIR"
done


