#!/bin/bash

# Define paths
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping
THREADS=16

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
    echo "  Mapping reads with minimap2..."
    minimap2 -ax sr -t $THREADS $REF_GENOME $R1 $R2 2> $LOG_FILE | \
    samtools view -@ $THREADS -b -o $INITIAL_BAM -
    
    # Check if mapping was successful
    if [ ! -s "$INITIAL_BAM" ]; then
        echo "  ERROR: Mapping failed for $sample"
        return 1
    fi
    
    # Step 2: Filter by mapping quality
    echo "  Filtering by mapping quality (MAPQ >= 20)..."
    samtools view -@ $THREADS -b -q 20 -o $FILTERED_BAM $INITIAL_BAM
    
    # Step 3: Sort by position
    echo "  Sorting by position..."
    samtools sort -@ $THREADS -o $SORTED_BAM $FILTERED_BAM
    
    # Step 4: Index the sorted BAM
    echo "  Indexing BAM file..."
    samtools index -@ $THREADS $SORTED_BAM
    
    # Step 5: Filter out potential PCR duplicates based on identical start-end positions
    echo "  Removing potential PCR duplicates..."
    samtools rmdup $SORTED_BAM $NODEDUP_BAM
    
    # Step 6: Index the deduplicated BAM file
    echo "  Indexing deduplicated BAM file..."
    samtools index -@ $THREADS $NODEDUP_BAM
    
    echo "  Processing complete for $sample"
}

# Process all samples in the reads directory
for R1_FILE in $READS_DIR/*_R1_minimal.fq.gz; do
    SAMPLE=$(basename $R1_FILE _R1_minimal.fq.gz)
    process_sample $SAMPLE
done

echo "Pipeline complete!"
