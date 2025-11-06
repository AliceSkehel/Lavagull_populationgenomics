#!/bin/bash

# Exit on error
set -e
set -o pipefail

REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna
R1=~/Sequencing_Combined/cleaned_fastq/LVGU_10_R1_cleaned.fq.gz
R2=~/Sequencing_Combined/cleaned_fastq/LVGU_10_R2_cleaned.fq.gz
OUT_DIR=~/Sequencing_Combined/complete_mapping/trial_and_error
SAMPLE=LVGU_10_testOct

# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

echo "Starting mapping pipeline for $SAMPLE..."
echo "Output directory: $OUT_DIR"

# Step 1: Map and filter
echo "[$(date)] Step 1: Mapping with minimap2..."
minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
  samtools view -bS -q 20 -@ 8 - > $OUT_DIR/${SAMPLE}_mapped.bam
echo "[$(date)] ✓ Mapped reads saved: ${SAMPLE}_mapped.bam"
ls -lh $OUT_DIR/${SAMPLE}_mapped.bam

# Step 2: Sort by name for fixmate
echo "[$(date)] Step 2: Sorting by name..."
samtools sort -n -@ 8 $OUT_DIR/${SAMPLE}_mapped.bam -o $OUT_DIR/${SAMPLE}_namesorted.bam
echo "[$(date)] ✓ Name-sorted BAM saved: ${SAMPLE}_namesorted.bam"
ls -lh $OUT_DIR/${SAMPLE}_namesorted.bam

# Step 3: Add mate tags
echo "[$(date)] Step 3: Running fixmate..."
samtools fixmate -m -@ 8 $OUT_DIR/${SAMPLE}_namesorted.bam $OUT_DIR/${SAMPLE}_fixmate.bam
echo "[$(date)] ✓ Fixmate BAM saved: ${SAMPLE}_fixmate.bam"
ls -lh $OUT_DIR/${SAMPLE}_fixmate.bam

# Step 4: Sort by coordinate
echo "[$(date)] Step 4: Sorting by coordinate..."
samtools sort -@ 8 $OUT_DIR/${SAMPLE}_fixmate.bam -o $OUT_DIR/${SAMPLE}_sorted.bam
echo "[$(date)] ✓ Coordinate-sorted BAM saved: ${SAMPLE}_sorted.bam"
ls -lh $OUT_DIR/${SAMPLE}_sorted.bam

# Step 5: Mark/remove duplicates
echo "[$(date)] Step 5: Marking duplicates..."
samtools markdup -r -@ 8 $OUT_DIR/${SAMPLE}_sorted.bam $OUT_DIR/${SAMPLE}.bam
echo "[$(date)] ✓ Final BAM saved: ${SAMPLE}.bam"
ls -lh $OUT_DIR/${SAMPLE}.bam

# Step 6: Index
echo "[$(date)] Step 6: Indexing..."
samtools index -@ 8 $OUT_DIR/${SAMPLE}.bam
echo "[$(date)] ✓ Index created: ${SAMPLE}.bam.bai"

# Optional: Clean up intermediate files
echo "[$(date)] Cleaning up intermediate files..."
read -p "Delete intermediate files? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm $OUT_DIR/${SAMPLE}_mapped.bam
    rm $OUT_DIR/${SAMPLE}_namesorted.bam
    rm $OUT_DIR/${SAMPLE}_fixmate.bam
    rm $OUT_DIR/${SAMPLE}_sorted.bam
    echo "[$(date)] ✓ Cleanup complete"
fi

echo "[$(date)] ✓ Pipeline complete!"
samtools flagstat $OUT_DIR/${SAMPLE}.bam
