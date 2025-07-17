#!/bin/bash

# Define paths
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping
THREADS=16

# Create output directories
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/flagstats
mkdir -p $OUTPUT_DIR/coverage_stats
mkdir -p $OUTPUT_DIR/mapping_logs

# Record start time
start_time=$(date +%s)
echo "Starting mapping and deduplication pipeline at $(date)"
echo "Reference genome: $REF_GENOME"
echo "Using $THREADS threads"
echo "----------------------------------------"

# Function to process a sample (mapping + proper deduplication)
process_sample() {
    local sample=$1
    echo "Processing sample: $sample"
    
    # Define paths for this sample
    R1=$READS_DIR/${sample}_R1_minimal.fq.gz
    R2=$READS_DIR/${sample}_R2_minimal.fq.gz
    
    # Define output files
    INITIAL_BAM=$OUTPUT_DIR/${sample}.initial.bam
    NAMESORT_BAM=$OUTPUT_DIR/${sample}.namesort.bam
    FIXMATE_BAM=$OUTPUT_DIR/${sample}.fixmate.bam
    POSITIONSORT_BAM=$OUTPUT_DIR/${sample}.sorted.bam
    DEDUP_BAM=$OUTPUT_DIR/${sample}.dedup.bam
    LOG_FILE=$OUTPUT_DIR/mapping_logs/${sample}.log
    
    # Check if input files exist
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Error: Input files for $sample not found. Skipping."
        return 1
    fi

    # Step 1: Run minimap2 mapping
    echo "  Mapping reads with minimap2..."
    minimap2 -ax sr -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
        $REF_GENOME $R1 $R2 2>> $LOG_FILE | \
        samtools view -@ $THREADS -bS - > $INITIAL_BAM
    
    # Step 2: Sort by name for fixmate
    echo "  Sorting by name for fixmate..."
    samtools sort -n -@ $THREADS -o $NAMESORT_BAM $INITIAL_BAM 2>> $LOG_FILE
    
    # Step 3: Run fixmate to add mate coordinates
    echo "  Adding mate coordinates with fixmate..."
    samtools fixmate -m -@ $THREADS $NAMESORT_BAM $FIXMATE_BAM 2>> $LOG_FILE
    
    # Step 4: Sort by coordinate
    echo "  Sorting by coordinate..."
    samtools sort -@ $THREADS -o $POSITIONSORT_BAM $FIXMATE_BAM 2>> $LOG_FILE
    
    # Step 5: Mark duplicates
    echo "  Marking duplicates..."
    samtools markdup -@ $THREADS $POSITIONSORT_BAM $DEDUP_BAM 2>> $LOG_FILE
    
    # Step 6: Index the final BAM
    echo "  Indexing BAM..."
    samtools index -@ $THREADS $DEDUP_BAM 2>> $LOG_FILE
    
    # Step 7: Generate mapping statistics
    echo "  Generating mapping statistics..."
    samtools flagstat -@ $THREADS $DEDUP_BAM > $OUTPUT_DIR/flagstats/${sample}.flagstat.txt
    
    # Step 8: Generate coverage statistics
    echo "  Calculating coverage statistics..."
    samtools coverage -@ $THREADS $DEDUP_BAM > $OUTPUT_DIR/coverage_stats/${sample}.coverage.txt
    
    # Calculate mean coverage
    mean_cov=$(samtools depth $DEDUP_BAM | awk '{sum+=$3} END {print sum/NR}')
    echo "  Mean coverage depth: $mean_cov"
    
    # Calculate duplication rate
    initial_reads=$(samtools view -c $INITIAL_BAM)
    final_reads=$(samtools view -c $DEDUP_BAM)
    dup_reads=$((initial_reads - final_reads))
    dup_percent=$(echo "scale=2; 100 * $dup_reads / $initial_reads" | bc)
    echo "  Duplication rate: ${dup_percent}% (${dup_reads} reads marked as duplicates)"
    
    # Clean up intermediate files
    echo "  Cleaning up intermediate files..."
    rm $INITIAL_BAM
    rm $NAMESORT_BAM
    rm $FIXMATE_BAM
    
    echo "  Completed processing of $sample"
    echo "----------------------------------------"
}

# Process all samples sequentially
echo "Detecting samples..."
samples=$(ls $READS_DIR/*_R1_minimal.fq.gz | sed 's/.*\///' | sed 's/_R1_minimal.fq.gz//')
echo "Found samples: $samples"
echo "----------------------------------------"

for sample in $samples; do
    process_sample $sample
done

# Create a summary report
echo "Creating summary report..."
echo "Mapping and Deduplication Summary" > $OUTPUT_DIR/pipeline_summary.txt
echo "===============================" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Date: $(date)" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Reference genome: $REF_GENOME" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Number of samples processed: $(echo $samples | wc -w)" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Threads used per sample: $THREADS" >> $OUTPUT_DIR/pipeline_summary.txt
echo "" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Mapping Statistics:" >> $OUTPUT_DIR/pipeline_summary.txt
echo "----------------" >> $OUTPUT_DIR/pipeline_summary.txt

for sample in $samples; do
    echo "Sample: $sample" >> $OUTPUT_DIR/pipeline_summary.txt
    
    # Extract key metrics
    if [ -f "$OUTPUT_DIR/flagstats/${sample}.flagstat.txt" ]; then
        total=$(grep "in total" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | cut -d' ' -f1)
        mapped=$(grep "mapped (" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | head -1 | cut -d' ' -f1)
        mapping_rate=$(grep "mapped (" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | head -1 | grep -o "[0-9.]*%" | head -1)
        properly_paired=$(grep "properly paired" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | cut -d' ' -f1)
        proper_rate=$(grep "properly paired" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | grep -o "[0-9.]*%" | head -1)
        dups=$(grep "duplicates" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | head -1 | cut -d' ' -f1)
        dup_rate=$(echo "scale=2; 100 * $dups / $total" | bc)
        
        echo "  Total reads: $total" >> $OUTPUT_DIR/pipeline_summary.txt
        echo "  Mapped reads: $mapped ($mapping_rate)" >> $OUTPUT_DIR/pipeline_summary.txt
        echo "  Properly paired: $properly_paired ($proper_rate)" >> $OUTPUT_DIR/pipeline_summary.txt
        echo "  Duplicates: $dups (${dup_rate}%)" >> $OUTPUT_DIR/pipeline_summary.txt
        
        # Add coverage information if available
        if [ -f "$OUTPUT_DIR/coverage_stats/${sample}.coverage.txt" ]; then
            mean_cov=$(awk 'NR>1 {sum+=$7; count++} END {print sum/count}' $OUTPUT_DIR/coverage_stats/${sample}.coverage.txt)
            cov_1x=$(awk 'NR>1 && $7>=1 {count++} END {print count/NR*100}' $OUTPUT_DIR/coverage_stats/${sample}.coverage.txt)
            cov_10x=$(awk 'NR>1 && $7>=10 {count++} END {print count/NR*100}' $OUTPUT_DIR/coverage_stats/${sample}.coverage.txt)
            
            echo "  Mean coverage: ${mean_cov}x" >> $OUTPUT_DIR/pipeline_summary.txt
            echo "  Genome covered at ≥1x: ${cov_1x}%" >> $OUTPUT_DIR/pipeline_summary.txt
            echo "  Genome covered at ≥10x: ${cov_10x}%" >> $OUTPUT_DIR/pipeline_summary.txt
        fi
    else
        echo "  No statistics available" >> $OUTPUT_DIR/pipeline_summary.txt
    fi
    
    echo "" >> $OUTPUT_DIR/pipeline_summary.txt
done

# Calculate and record runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))

echo "Runtime: ${hours}h ${minutes}m ${seconds}s" >> $OUTPUT_DIR/pipeline_summary.txt
echo "" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Final BAM files are in: $OUTPUT_DIR" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Flagstat files are in: $OUTPUT_DIR/flagstats" >> $OUTPUT_DIR/pipeline_summary.txt
echo "Coverage statistics are in: $OUTPUT_DIR/coverage_stats" >> $OUTPUT_DIR/pipeline_summary.txt

echo "Pipeline complete!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "Summary report is available at: $OUTPUT_DIR/pipeline_summary.txt"
