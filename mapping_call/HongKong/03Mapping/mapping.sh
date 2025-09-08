#!/bin/bash
set -eu  

# Path to the reference genome file
reference=../../../data/Reference/data/GCF_015237465.2/Chelonia_mydas_GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

# Directory containing FASTQ files
PACIFIC_DATA_DIR=../02Ready_to_mapping/fastqgz


# Output directories
SAM_DIR="sam"
BWA_DIR="bwa"
STATS_DIR="stats"

# Number of CPU cores to use
CORE=5

# Create output directories if they do not exist
mkdir -p "$SAM_DIR" "$BWA_DIR" "$STATS_DIR"


# Loop through each pair of FASTQ files for mapping and processing
for fastq_file_1 in "$PACIFIC_DATA_DIR"/*_1.fastq.gz; do
    # Extract the sample name from the file name
    sample_name=$(basename "$fastq_file_1" _1.fastq.gz)
    
    # Define the reverse read file
    fastq_file_2="$PACIFIC_DATA_DIR/${sample_name}_2.fastq.gz"
    
    # Check if the reverse read file exists
    if [[ ! -f "$fastq_file_2" ]]; then
        echo "Warning: Paired file for $fastq_file_1 not found. Skipping..."
        continue
    fi
    
    # Mapping reads to the reference genome with BWA
    echo "Mapping $sample_name to the reference..."
    bwa mem -t $CORE -R "@RG\tID:$sample_name\tSM:$sample_name\tPL:illumina\tLB:${sample_name}_library\tPU:unit1" \
    "$reference" "$fastq_file_1" "$fastq_file_2" > "$SAM_DIR/$sample_name.sam"
    
    # Convert SAM to sorted BAM
    echo "Converting $sample_name.sam to $sample_name.bam..."
    samtools sort -@ $CORE "$SAM_DIR/$sample_name.sam" -o "$BWA_DIR/$sample_name.bam"

    # Index the BAM file
    echo "Indexing $sample_name.bam..."
    samtools index "$BWA_DIR/$sample_name.bam"

    # Generate alignment statistics
    echo "Generating alignment stats for $sample_name..."
    samtools stats -@ $CORE "$BWA_DIR/$sample_name.bam" > "$STATS_DIR/$sample_name.stats"
done

echo "Mapping, conversion, and indexing completed."
