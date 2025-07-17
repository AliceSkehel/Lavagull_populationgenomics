#!/bin/bash

conda activate LVGU1

# Define paths
FASTQ_DIR=~/Sequencing_Round1
FASTQC_RESULTS_DIR=~/Sequencing_Round1/fastqc_results

# Navigate to the directory with FASTQ files
cd $FASTQ_DIR

# Run FastQC on all FASTQ files
echo "Running FastQC on all files..."
fastqc -t 8 -o $FASTQC_RESULTS_DIR *.fq.gz

# Run a multiqc on all the FASTQC files to generate just one report
multiqc fastqc_results -o multiqc_report

# Copy results to the local machine by going into a terminal that's not logged onto the server and typing this in:
# scp -r askehel@alice:~/Sequencing_Round1/multiqc_report/* /Users/aliceskehel/VSCode/MultiQC_Round1_Report/
