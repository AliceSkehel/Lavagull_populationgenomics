#!/bin/bash

# Mapping Script from the 9th June
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping

mkdir -p $OUTPUT_DIR

for R1_FILE in $READS_DIR/*_R1_minimal.fq.gz; do
    SAMPLE=$(basename $R1_FILE _R1_minimal.fq.gz)
    R1=$READS_DIR/${SAMPLE}_R1_minimal.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_minimal.fq.gz
    
    echo "Processing: $SAMPLE"

  # Skip if already processed
    if [ -f "$OUTPUT_DIR/${SAMPLE}.bam" ] && [ -f "$OUTPUT_DIR/${SAMPLE}.bam.bai" ]; then
        echo "  $SAMPLE already completed, skipping..."
        continue
    fi
    
    (
    # Map with minimap2, filter by quality, sort by name for fixmate
    minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
    samtools view -bS -q 20 -@ 8 - | \
    samtools sort -n -@ 8 - | \
    samtools fixmate -m -@ 8 - - | \
    samtools sort -@ 8 - | \
samtools view -h -@ 8 - | \
    awk '/^@/ || /MC:Z:/' | \
    samtools view -b -@ 8 - | \
    samtools markdup -r -@ 8 - $OUTPUT_DIR/${SAMPLE}.bam
    
    samtools index -@ 8 $OUTPUT_DIR/${SAMPLE}.bam
    ) &
done

wait

    echo "  ✅ Completed: ${SAMPLE}.bam (deduplicated and indexed)"
    ) &
done

wait

echo "🎉 All samples processed successfully!"
echo ""
echo "Output files in: $OUTPUT_DIR"
echo "Files are ready for variant calling, Qualimap, and other analyses."
