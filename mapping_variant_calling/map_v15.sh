#!/bin/bash

# Simple Mapping Script
READS_DIR=~/Sequencing_Round2/fastp_polyg_cleaned/cleaned_fastq
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.fna
OUTPUT_DIR=~/Sequencing_Round2/complete_mapping/unfiltered_ref

mkdir -p $OUTPUT_DIR

for R1_FILE in $READS_DIR/*_R1_cleaned.fq.gz; do
    SAMPLE=$(basename $R1_FILE _R1_cleaned.fq.gz)
    
    echo "Processing: $SAMPLE"
    
    # Skip if already done
    if [ -f "$OUTPUT_DIR/${SAMPLE}.bam" ]; then
        echo "  Already done, skipping..."
        continue
    fi
    
    R1=$READS_DIR/${SAMPLE}_R1_cleaned.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_cleaned.fq.gz
    
    # Map and process
    minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
    samtools view -bS -q 20 -@ 8 - | \
    samtools sort -n -@ 8 - | \
    samtools fixmate -m -@ 8 - - | \
    samtools sort -@ 8 - | \
    samtools view -h -@ 8 - | \
    awk '/^@/ || /MC:Z:/' | \
    samtools view -b -@ 8 - | \
    samtools markdup -r -@ 8 - $OUTPUT_DIR/${SAMPLE}.bam
    
    # Index
    samtools index -@ 8 $OUTPUT_DIR/${SAMPLE}.bam
    
    echo "  Done: ${SAMPLE}.bam"
done

echo "All samples processed!"
