#!/bin/bash

# Updated Mapping Script for Round2 Poly-G Cleaned Data
READS_DIR=~/Sequencing_Round2/fastp_polyg_cleaned/cleaned_fastq
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round2/complete_mapping

mkdir -p $OUTPUT_DIR

echo "====================================="
echo "Mapping Pipeline - Round 2 Data"
echo "Input: $READS_DIR"
echo "Reference: $REF_GENOME"
echo "Output: $OUTPUT_DIR"
echo "====================================="

for R1_FILE in $READS_DIR/*_R1_cleaned.fq.gz; do
    # Extract sample name (LVGU_X from LVGU_X_R1_cleaned.fq.gz)
    SAMPLE=$(basename $R1_FILE _R1_cleaned.fq.gz)
    R1=$READS_DIR/${SAMPLE}_R1_cleaned.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_cleaned.fq.gz
    
    echo "Processing: $SAMPLE"

    # Skip if already processed
    if [ -f "$OUTPUT_DIR/${SAMPLE}.bam" ] && [ -f "$OUTPUT_DIR/${SAMPLE}.bam.bai" ]; then
        echo "  $SAMPLE already completed, skipping..."
        continue
    fi
    
    # Check if input files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "  ❌ Missing input files for $SAMPLE"
        continue
    fi
    
    echo "  Input files:"
    echo "    R1: $R1"
    echo "    R2: $R2"
    
    (
    # Map with minimap2, filter by quality, sort by name for fixmate
    echo "  Starting mapping pipeline for $SAMPLE..."
    minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
    samtools view -bS -q 20 -@ 8 - | \
    samtools sort -n -@ 8 - | \
    samtools fixmate -m -@ 8 - - | \
    samtools sort -@ 8 - | \
    samtools view -h -@ 8 - | \
    awk '/^@/ || /MC:Z:/' | \
    samtools view -b -@ 8 - | \
    samtools markdup -r -@ 8 - $OUTPUT_DIR/${SAMPLE}.bam
    
    # Index the BAM file
    samtools index -@ 8 $OUTPUT_DIR/${SAMPLE}.bam
    
    echo "  ✅ Completed: ${SAMPLE}.bam (deduplicated and indexed)"
    ) &
done

wait

echo ""
echo "🎉 All samples processed successfully!"
echo ""
echo "Output files in: $OUTPUT_DIR"
echo "Files are ready for variant calling, Qualimap, and other analyses."
echo ""
echo "Generated files:"
ls -lh $OUTPUT_DIR/*.bam 2>/dev/null | head -5
echo ""
