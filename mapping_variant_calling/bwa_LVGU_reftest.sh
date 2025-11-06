#!/bin/bash

# BWA mapping to LVGU_60 - using map_v16 pipeline
# Test samples: LVGU_18, LVGU_7, LVGU_15, LVGU_56, LVGU_10

READS_DIR=~/Sequencing_Combined/cleaned_fastq
REF_GENOME=~/Sequencing_Mapping/LVGU_60_ref.fasta
OUTPUT_DIR=~/Sequencing_Mapping/outputs/bwa_bams

mkdir -p $OUTPUT_DIR

# Index reference
bwa index $REF_GENOME

# Process each test sample
for SAMPLE in LVGU_18 LVGU_7 LVGU_15 LVGU_56 LVGU_10; do
    echo "Processing: $SAMPLE"
    
    if [ -f "$OUTPUT_DIR/${SAMPLE}.bam" ]; then
        echo "  Already done, skipping..."
        continue
    fi
    
    R1=$READS_DIR/${SAMPLE}_R1_cleaned.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_cleaned.fq.gz
    
    bwa mem -t 8 $REF_GENOME $R1 $R2 | \
    samtools view -bS -q 20 -@ 8 - | \
    samtools sort -n -@ 8 - | \
    samtools fixmate -m -@ 8 - - | \
    samtools sort -@ 8 - | \
    samtools view -h -@ 8 - | \
    awk '/^@/ || /MC:Z:/' | \
    samtools view -b -@ 8 - | \
    samtools markdup -r -@ 8 - $OUTPUT_DIR/${SAMPLE}.bam
    
    samtools index -@ 8 $OUTPUT_DIR/${SAMPLE}.bam
    
    echo "  Done: ${SAMPLE}.bam"
done

echo "All samples processed!"
