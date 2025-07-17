#!/bin/bash

# Define paths
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping
THREADS=4  # Reduce threads per sample since we're running in parallel

# Maximum number of parallel jobs to run
MAX_PARALLEL_JOBS=4  # Adjust based on your system's capabilities

# Create output directories
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/mapping_logs
mkdir -p $OUTPUT_DIR/dedup

# Function to process a sample (mapping only)
process_sample() {
    local sample=$1
    echo "Processing sample: $sample"
    
    # Define paths for this sample
    R1=$READS_DIR/${sample}_R1_minimal.fq.gz
    R2=$READS_DIR/${sample}_R2_minimal.fq.gz
    
    # Define output files
    INITIAL_BAM=$OUTPUT_DIR/${sample}.initial.bam
    FILTERED_BAM=$OUTPUT_DIR/${sample}.filtered.bam
    SORTED_BAM=$OUTPUT_DIR/${sample}.sorted.bam
    NODEDUP_BAM=$OUTPUT_DIR/dedup/${sample}.nodedup.bam
    LOG_FILE=$OUTPUT_DIR/mapping_logs/${sample}.log
    
    # Step 1: Run minimap2 mapping
    echo "  Mapping reads with minimap2 for $sample..."
    minimap2 -ax sr -t $THREADS $REF_GENOME $R1 $R2 2> $LOG_FILE | \
    samtools view -@ $THREADS -b -o $INITIAL_BAM -
    
    # Check if mapping was successful
    if [ ! -s "$INITIAL_BAM" ]; then
        echo "  ERROR: Mapping failed for $sample"
        return 1
    fi
    
    # Step 2: Filter by mapping quality
    echo "  Filtering by mapping quality (MAPQ >= 20) for $sample..."
    samtools view -@ $THREADS -b -q 20 -o $FILTERED_BAM $INITIAL_BAM
    
    # Step 3: Sort by position
    echo "  Sorting by position for $sample..."
    samtools sort -@ $THREADS -o $SORTED_BAM $FILTERED_BAM
    
    # Step 4: Index the sorted BAM
    echo "  Indexing BAM file for $sample..."
    samtools index -@ $THREADS $SORTED_BAM
    
    # Step 5: Filter out potential PCR duplicates based on identical start-end positions
    echo "  Removing potential PCR duplicates for $sample..."
    samtools rmdup $SORTED_BAM $NODEDUP_BAM
    
    # Step 6: Index the deduplicated BAM file
    echo "  Indexing deduplicated BAM file for $sample..."
    samtools index -@ $THREADS $NODEDUP_BAM
    
    echo "  Processing complete for $sample"
}

# Function to run jobs in parallel
run_parallel() {
    # Create an array to store process IDs
    local pids=()
    
    # Process all samples in the reads directory
    for R1_FILE in $READS_DIR/*_R1_minimal.fq.gz; do
        local SAMPLE=$(basename $R1_FILE _R1_minimal.fq.gz)
        
        # Skip if output already exists (optional)
        if [ -f "$OUTPUT_DIR/dedup/${SAMPLE}.nodedup.bam" ]; then
            echo "Sample $SAMPLE already processed, skipping..."
            continue
        fi
        
        # Run process_sample in the background
        process_sample $SAMPLE &
        
        # Store the process ID
        pids+=($!)
        
        # Check if we've reached the maximum number of parallel jobs
        if [ ${#pids[@]} -ge $MAX_PARALLEL_JOBS ]; then
            # Wait for one of the jobs to finish
            wait -n
            
            # Remove completed jobs from the array
            local new_pids=()
            for pid in "${pids[@]}"; do
                if kill -0 $pid 2>/dev/null; then
                    new_pids+=($pid)
                fi
            done
            pids=("${new_pids[@]}")
        fi
    done
    
    # Wait for all remaining jobs to complete
    wait
}

echo "Starting parallel processing with $MAX_PARALLEL_JOBS simultaneous jobs..."
run_parallel
echo "Pipeline complete!"
