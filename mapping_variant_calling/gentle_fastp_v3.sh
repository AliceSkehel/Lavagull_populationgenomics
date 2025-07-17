#!/bin/bash
# Define paths with distinct names for this gentler cleaning
INPUT_DIR=~/Sequencing_Round1
OUTPUT_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REPORT_DIR=~/Sequencing_Round1/fastp_minimal_reports
# Create output directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $REPORT_DIR
# Process all samples using a more dynamic approach
cd $INPUT_DIR
# Find all unique sample prefixes
echo "Detecting samples..."
samples=$(ls *_L1_1.fq.gz | sed 's/_CKDL.*$//')
echo "Found samples: $samples"
# Process each sample
for sample in $samples; do
    echo "======================================"
    echo "Processing sample: $sample"
    echo "======================================"
    
    r1="${sample}_CKDL250009306-1A_22THFLLT4_L1_1.fq.gz"
    r2="${sample}_CKDL250009306-1A_22THFLLT4_L1_2.fq.gz"
    
    # Verify files exist
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Error: Input files for $sample not found. Skipping."
        continue
    fi
    
    # Run fastp with MINIMAL parameters spread across multiple lines
    # IMPORTANT: No spaces or comments after the backslashes!
    fastp \
-i "$r1" \
-I "$r2" \
-o "$OUTPUT_DIR/${sample}_R1_minimal.fq.gz" \
-O "$OUTPUT_DIR/${sample}_R2_minimal.fq.gz" \
-q 20 \
--adapter_sequence=GATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--detect_adapter_for_pe \
--poly_g_min_len 10 \
--thread 8 \
-h "$REPORT_DIR/${sample}_fastp_minimal.html" \
-j "$REPORT_DIR/${sample}_fastp_minimal.json" \
--verbose
    
    # Check if fastp completed successfully
    if [ $? -eq 0 ]; then
        echo "Successfully processed $sample"
    else
        echo "Error processing $sample"
        echo "Continuing with next sample..."
    fi
done
echo "All samples processed!"
# Run FastQC on minimally cleaned files to verify quality
echo "Running FastQC on minimally cleaned files..."
mkdir -p $OUTPUT_DIR/fastqc_minimal_results
fastqc -t 8 -o $OUTPUT_DIR/fastqc_minimal_results $OUTPUT_DIR/*.fq.gz
# Run MultiQC on FastQC results for easy comparison
echo "Running MultiQC..."
cd $OUTPUT_DIR
multiqc fastqc_minimal_results -o multiqc_minimal_report
echo "====================================="
echo "Processing Summary:"
echo "Total samples processed: $(echo $samples | wc -w)"
echo "Minimally cleaned files are in: $OUTPUT_DIR"
echo "FASTP reports are in: $REPORT_DIR"
echo "FastQC results are in: $OUTPUT_DIR/fastqc_minimal_results"
echo "MultiQC report is in: $OUTPUT_DIR/multiqc_minimal_report"
echo "====================================="
