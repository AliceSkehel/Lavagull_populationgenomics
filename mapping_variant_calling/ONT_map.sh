#!/bin/bash

# Simple ONT mtDNA Mapping Script
READS_DIR=~/mitogenome/mito_fastq_merged
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.fasta
OUTPUT_DIR=~/mitogenome/mapped_bams

mkdir -p $OUTPUT_DIR

echo "Starting ONT mtDNA mapping..."

for FASTQ_FILE in $READS_DIR/LVGU_*.fastq.gz; do
    SAMPLE=$(basename $FASTQ_FILE .fastq.gz)
    
    echo "Processing: $SAMPLE"
    
    # Map and sort
    minimap2 -ax map-ont -t 8 $REF_GENOME $FASTQ_FILE | \
    samtools sort -@ 8 -o $OUTPUT_DIR/${SAMPLE}_ONT.bam
    
    # Index
    samtools index $OUTPUT_DIR/${SAMPLE}_ONT.bam
    
    echo "  ✅ Done: ${SAMPLE}_ONT.bam"
done

echo "🎉 All samples completed!"
ls -lh $OUTPUT_DIR/*_ONT.bam
