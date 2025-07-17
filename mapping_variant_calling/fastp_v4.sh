#!/bin/bash

# Poly-G and Adapter Cleanup on Already Filtered Files

# Define paths
INPUT_DIR=~/Sequencing_Round2/fastp_processed/filtered_fastq
OUTPUT_DIR=~/Sequencing_Round2/fastp_polyg_cleaned
CLEANED_DIR=$OUTPUT_DIR/cleaned_fastq

# Create output directories
mkdir -p $CLEANED_DIR

echo "====================================="
echo "Poly-G and Adapter Cleanup Pipeline"
echo "Input (already filtered): $INPUT_DIR"
echo "Output (poly-G cleaned): $CLEANED_DIR"
echo "====================================="

cd $INPUT_DIR

# Count input files
TOTAL_SAMPLES=$(ls *_R1_filtered.fq.gz 2>/dev/null | wc -l)
echo "Found $TOTAL_SAMPLES samples to clean..."
echo ""

# Process each sample
for r1_file in *_R1_filtered.fq.gz; do
    if [ -f "$r1_file" ]; then
        # Extract sample name (remove _R1_filtered.fq.gz)
        sample_name=$(basename "$r1_file" _R1_filtered.fq.gz)
        r2_file="${sample_name}_R2_filtered.fq.gz"
        
        echo "Processing $sample_name..."
        
        if [ -f "$r2_file" ]; then
            # Run FASTP with ONLY adapter removal and poly-G trimming
            fastp \
                -i "$r1_file" \
                -I "$r2_file" \
                -o "$CLEANED_DIR/${sample_name}_R1_cleaned.fq.gz" \
                -O "$CLEANED_DIR/${sample_name}_R2_cleaned.fq.gz" \
                --detect_adapter_for_pe \
                --trim_poly_g \
                --poly_g_min_len 10 \
                --disable_quality_filtering \
                --disable_length_filtering
            
            echo "✅ $sample_name completed"
        else
            echo "❌ R2 file not found for $sample_name"
        fi
    fi
done

echo ""
echo "🎉 Poly-G cleanup completed!"
echo "Cleaned files: $CLEANED_DIR"
echo ""

# Run FastQC and MultiQC on cleaned files
echo "====================================="
echo "Running FastQC on poly-G cleaned files..."
echo "====================================="

# Define FastQC output directory
FASTQC_OUTPUT=$OUTPUT_DIR/fastqc_cleaned_reports
MULTIQC_OUTPUT=$OUTPUT_DIR/multiqc_cleaned_summary

# Create FastQC output directory
mkdir -p $FASTQC_OUTPUT
mkdir -p $MULTIQC_OUTPUT

# Count cleaned files
TOTAL_FILES=$(ls $CLEANED_DIR/*.fq.gz 2>/dev/null | wc -l)
echo "Found $TOTAL_FILES cleaned FASTQ files to analyze..."
echo ""

# Run FastQC on all cleaned files
fastqc $CLEANED_DIR/*.fq.gz -o $FASTQC_OUTPUT -t 4

echo ""
echo "FastQC completed! Generated $(ls $FASTQC_OUTPUT/*.html | wc -l) HTML reports"
echo ""

# Generate MultiQC summary report
echo "Creating MultiQC summary report..."
cd $FASTQC_OUTPUT
multiqc . --title "Poly_G_Cleaned_Data_QC_Summary" -o $MULTIQC_OUTPUT

echo ""
echo "====================================="
echo "🎉 Complete Cleanup Pipeline Finished!"
echo ""
echo "Results:"
echo "  Original filtered files: $INPUT_DIR"
echo "  Poly-G cleaned files: $CLEANED_DIR"
echo "  FastQC reports: $FASTQC_OUTPUT"
echo "  MultiQC summary: $MULTIQC_OUTPUT"
echo ""
echo "To view MultiQC report:"
echo "  firefox $MULTIQC_OUTPUT/multiqc_report.html"
echo ""
echo "File naming:"
echo "  Input: LVGU_X_R1_filtered.fq.gz"
echo "  Output: LVGU_X_R1_cleaned.fq.gz"
echo "====================================="
