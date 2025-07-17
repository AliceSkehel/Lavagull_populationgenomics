#!/bin/bash

# Streamlined mapping pipeline - no intermediate unsorted files!

# Define paths
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping
THREADS=4
MAX_PARALLEL_JOBS=4

# Create output directories
mkdir -p $OUTPUT_DIR/logs

# Function to process a sample
process_sample() {
    local sample=$1
    echo "Processing sample: $sample"
    
    # Define file paths
    R1=$READS_DIR/${sample}_R1_minimal.fq.gz
    R2=$READS_DIR/${sample}_R2_minimal.fq.gz
    FINAL_BAM=$OUTPUT_DIR/${sample}.final.bam
    LOG_FILE=$OUTPUT_DIR/logs/${sample}.log
    
    # Skip if already processed
    if [ -f "$FINAL_BAM.bai" ]; then
        echo "  Sample $sample already processed, skipping..."
        return 0
    fi
    
    echo "  Mapping, filtering, sorting, and deduplicating $sample..."
    
    # ONE-STEP PIPELINE: Map → Filter → Sort → Deduplicate → Index
    minimap2 -ax sr -t $THREADS $REF_GENOME $R1 $R2 2> $LOG_FILE | \
    samtools view -@ $THREADS -b -q 20 - | \
    samtools sort -@ $THREADS - | \
    samtools markdup -r -@ $THREADS - $FINAL_BAM
    
    # Check if pipeline succeeded
    if [ $? -ne 0 ] || [ ! -s "$FINAL_BAM" ]; then
        echo "  ERROR: Pipeline failed for $sample"
        rm -f $FINAL_BAM  # Clean up partial file
        return 1
    fi
    
    # Index the final BAM
    echo "  Indexing $sample..."
    samtools index -@ $THREADS $FINAL_BAM
    
    echo "  ✅ Processing complete for $sample"
}

# Function to run jobs in parallel
run_parallel() {
    local pids=()
    
    for R1_FILE in $READS_DIR/*_R1_minimal.fq.gz; do
        local SAMPLE=$(basename $R1_FILE _R1_minimal.fq.gz)
        
        # Run in background
        process_sample $SAMPLE &
        pids+=($!)
        
        # Manage parallel job limit
        if [ ${#pids[@]} -ge $MAX_PARALLEL_JOBS ]; then
            wait -n  # Wait for one job to finish
            # Clean up completed jobs from array
            local new_pids=()
            for pid in "${pids[@]}"; do
                if kill -0 $pid 2>/dev/null; then
                    new_pids+=($pid)
                fi
            done
            pids=("${new_pids[@]}")
        fi
    done
    
    # Wait for remaining jobs
    wait
}

echo "🚀 Starting streamlined pipeline with $MAX_PARALLEL_JOBS parallel jobs..."
echo "📁 Output: Ready-to-use BAM files will be saved as *.final.bam"
run_parallel
echo "🎉 Pipeline complete! All files are sorted, deduplicated, and indexed."
