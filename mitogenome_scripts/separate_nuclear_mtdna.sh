#!/bin/bash

# Extract mtDNA
INPUT_DIR=~/Sequencing_Round2/complete_mapping/subset
OUTPUT_DIR=~/Sequencing_Round2/mtdna_only

mkdir -p $OUTPUT_DIR

for BAM_FILE in $INPUT_DIR/*.bam; do
    SAMPLE=$(basename $BAM_FILE .bam)
    echo "Processing: $SAMPLE"
    
    samtools view -b $BAM_FILE OZ118781.1 > $OUTPUT_DIR/${SAMPLE}_mtdna.bam
    samtools index $OUTPUT_DIR/${SAMPLE}_mtdna.bam
done

echo "Done!"
