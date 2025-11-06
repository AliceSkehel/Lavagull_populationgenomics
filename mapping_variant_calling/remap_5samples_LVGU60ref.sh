#!/bin/bash

# Remap 5 samples to LVGU_60 reference using working minimap2 pipeline
# Using brian123 environment with samtools 1.19

REF=~/Sequencing_Mapping/LVGU_60_ref.fasta
READS_DIR=~/Sequencing_Combined/cleaned_fastq
OUT_DIR=~/Sequencing_Combined/lavagull_ref_minimap2
mkdir -p $OUT_DIR

SAMPLES=(LVGU_5 LVGU_7 LVGU_13 LVGU_15 LVGU_18 LVGU_37 LVGU_42)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."
    
    R1=$READS_DIR/${SAMPLE}_R1_cleaned.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_cleaned.fq.gz
    
    minimap2 -ax sr -t 8 $REF $R1 $R2 | \
      samtools view -bS -q 20 -@ 8 - | \
      samtools sort -n -@ 8 - | \
      samtools fixmate -m -@ 8 - - | \
      samtools sort -@ 8 - | \
      samtools markdup -r -@ 8 - $OUT_DIR/${SAMPLE}.LavaGull_minimap2.bam
    
    samtools index $OUT_DIR/${SAMPLE}.LavaGull_minimap2.bam
    
    echo "Done: $SAMPLE"
done

echo "All 7 samples remapped to LVGU_60 reference"
