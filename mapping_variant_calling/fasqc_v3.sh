#!/bin/bash

# Define paths
FILTERED_DIR=~/Sequencing_Round2/fastp_processed/filtered_fastq
FASTQC_OUTPUT=~/Sequencing_Round2/fastqc_filtered_reports
MULTIQC_OUTPUT=~/Sequencing_Round2/multiqc_filtered_summary

# Create output directories
echo "Creating output directories..."
mkdir -p $FASTQC_OUTPUT
mkdir -p $MULTIQC_OUTPUT

echo "====================================="
echo "FastQC Analysis on Filtered Files"
echo "Input: $FILTERED_DIR"
echo "FastQC output: $FASTQC_OUTPUT"
echo "MultiQC output: $MULTIQC_OUTPUT"
echo "====================================="

# Count files
TOTAL_FILES=$(ls $FILTERED_DIR/*.fq.gz 2>/dev/null | wc -l)
echo "Found $TOTAL_FILES filtered FASTQ files to analyze..."
echo ""

# Run FastQC on all filtered files
echo "Running FastQC on all filtered files..."
fastqc $FILTERED_DIR/*.fq.gz -o $FASTQC_OUTPUT -t 4

echo ""
echo "FastQC completed! Generated $(ls $FASTQC_OUTPUT/*.html | wc -l) HTML reports"
echo ""

# Generate MultiQC summary report
echo "Creating MultiQC summary report..."
cd $FASTQC_OUTPUT
multiqc . --title "Filtered_Data_QC_Summary" -o $MULTIQC_OUTPUT

echo ""
echo "====================================="
echo "🎉 Quality Control Analysis Complete!"
echo ""
echo "Results:"
echo "  FastQC reports: $FASTQC_OUTPUT/"
echo "  MultiQC summary: $MULTIQC_OUTPUT/"
echo ""
echo "To view results:"
echo "  firefox $MULTIQC_OUTPUT/multiqc_report.html"
echo "====================================="
