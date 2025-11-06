#!/bin/bash

# Gentle fastp filtering for seq_batch1 and seq_batch2
# Using less stringent parameters to preserve more reads

INPUT_DIR="/home/askehel/Sequencing_streamlined/filtering_fastp"
OUTPUT_DIR="/home/askehel/Sequencing_streamlined/fastp_cleaned"
REPORT_DIR="/home/askehel/Sequencing_streamlined/fastp_reports"

mkdir -p $OUTPUT_DIR $REPORT_DIR

# Process all samples in seq_batch1 and seq_batch2
for batch in seq_batch1 seq_batch2; do
    for sample_dir in $INPUT_DIR/$batch/LVGU_*; do
        sample=$(basename $sample_dir)
        
        # Find R1 and R2 files
        r1=$(find $sample_dir -name "*_R1*" -o -name "*_1.fq.gz" | head -n 1)
        r2=$(find $sample_dir -name "*_R2*" -o -name "*_2.fq.gz" | head -n 1)
        
        echo "Processing $sample..."
        
        mkdir -p $OUTPUT_DIR/$sample
        
        fastp \
            -i $r1 \
            -I $r2 \
            -o $OUTPUT_DIR/$sample/${sample}_R1_cleaned.fq.gz \
            -O $OUTPUT_DIR/$sample/${sample}_R2_cleaned.fq.gz \
            -q 20 \
            --adapter_sequence=GATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --detect_adapter_for_pe \
            --poly_g_min_len 10 \
            --thread 12 \
            -h $REPORT_DIR/${sample}_report.html \
            -j $REPORT_DIR/${sample}_report.json
    done
done

echo "Done!"
echo "Cleaned files: $OUTPUT_DIR"
echo "Reports: $REPORT_DIR"
