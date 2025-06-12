#!/bin/bash

# レファレンスファイルのパスを指定
REFERENCE="/home/ofuji/Lab/seaturtle/reference/ncbi_dataset/data/GCF_015237465.2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna"

# fastq ファイルが保存されているディレクトリを指定
FASTQ_DIR="fastq"

# マッピング結果を保存するディレクトリを指定
MAPPING_DIR="sam"

# BWA ファイルを保存するディレクトリを指定
BWA_DIR="bwa"

# BWA インデックスファイルを保存するディレクトリを指定
BWA_INDEX_DIR="bwa_bai"

# STATS ファイルを保存するディレクトリを指定
STATS_DIR="stats"

# ディレクトリが存在しない場合は作成
mkdir -p "$MAPPING_DIR"
mkdir -p "$BWA_DIR"
mkdir -p "$BWA_INDEX_DIR"
mkdir -p "$STATS_DIR"

# fastq ファイルごとにmapping、mapping情報の保存、インデックス付けを行う
for fastq_file_1 in "$FASTQ_DIR"/*_1.fastq; do
    # サンプル名を取得
    sample_name=$(basename "$fastq_file_1" _1.fastq)
    
    # Reverseのfastqファイルを指定
    fastq_file_2="$FASTQ_DIR/${sample_name}_2.fastq"
    
    # マッピング
    echo "Mapping $sample_name to the reference..."
    bwa mem -t 2 "$REFERENCE" "$fastq_file_1" "$fastq_file_2" > "$MAPPING_DIR/$sample_name.sam"
    
    # SAM ファイルを BAM ファイルに変換
    echo "Converting $sample_name.sam to $sample_name.bam..."
    samtools sort -@ 2 "$MAPPING_DIR/$sample_name.sam" > "$MAPPING_DIR/$sample_name.bam"

    # BAM ファイルにインデックスを付ける
    echo "Indexing $sample_name.bam..."
    samtools index "$MAPPING_DIR/$sample_name.bam" > "$MAPPING_DIR/$sample_name.bam.bai"

    # STATS ファイルを作成
    echo "Generating alignment stats for $sample_name..."
    samtools stats -@ 2 "$MAPPING_DIR/$sample_name.bam" > "$STATS_DIR/$sample_name.stats"

    # BWA ファイルを bwa ディレクトリに移動
    mv "$MAPPING_DIR/$sample_name.bam" "$BWA_DIR/"

    # .bai ファイルを bwa_bai ディレクトリに移動
    mv "$MAPPING_DIR/$sample_name.bam.bai" "$BWA_INDEX_DIR/"
done

echo "Mapping, conversion, and indexing completed."
