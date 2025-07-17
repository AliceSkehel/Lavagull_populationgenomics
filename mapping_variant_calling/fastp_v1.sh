#!/bin/bash

# Simple FASTP Quality Filtering
conda activate LVGU1

# Define paths
RAW_FASTQ_DIR=Sequencing_Round2/fastq_Illumina_June/01.RawData
OUTPUT_DIR=Sequencing_Round2/fastp_processed
FILTERED_DIR=$OUTPUT_DIR/filtered_fastq

# Create output directories
mkdir -p $FILTERED_DIR

echo "Processing samples with FASTP..."
echo "Raw data: $RAW_FASTQ_DIR"
echo "Output: $FILTERED_DIR"

cd $RAW_FASTQ_DIR

# Process each sample
for sample_dir in LVGU_*/; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        echo "Processing $sample_name..."
        
        # Find R1 and R2 files (Novogene naming convention)
        R1=$(ls "$sample_dir"*_1.fq.gz 2>/dev/null | head -1)
        R2=$(ls "$sample_dir"*_2.fq.gz 2>/dev/null | head -1)
        
        if [ -f "$R1" ] && [ -f "$R2" ]; then
            # Simple FASTP with essential filters
            fastp \
                -i "$R1" \
                -I "$R2" \
                -o "$FILTERED_DIR/${sample_name}_R1_filtered.fq.gz" \
                -O "$FILTERED_DIR/${sample_name}_R2_filtered.fq.gz" \
                --qualified_quality_phred 20 \
                --cut_tail \
                --cut_tail_window_size 4 \
                --cut_tail_mean_quality 20 \
                --detect_adapter_for_pe
            
            echo "✅ $sample_name completed"
        else
            echo "❌ Files not found for $sample_name"
        fi
    fi
done

echo "🎉 FASTP processing completed!"
echo "Filtered files: $FILTERED_DIR"
