#!/bin/bash

# Define paths
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.fna
OUTPUT_DIR=~/Sequencing_Round1/minimap2_mappings
THREADS=16  # Increased to 16 threads

# Create output directory
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/flagstats
mkdir -p $OUTPUT_DIR/coverage_stats

# Record start time
start_time=$(date +%s)
echo "Starting minimap2 mapping at $(date)"
echo "Reference genome: $REF_GENOME"
echo "Using $THREADS threads"
echo "----------------------------------------"

# Function to process a sample
process_sample() {
    local sample=$1
    echo "Processing sample: $sample"
    
    # Define paths for this sample
    R1=$READS_DIR/${sample}_R1_minimal.fq.gz
    R2=$READS_DIR/${sample}_R2_minimal.fq.gz
    OUT_BAM=$OUTPUT_DIR/${sample}.sorted.bam
    
    # Check if input files exist
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Error: Input files for $sample not found. Skipping."
        return 1
    fi

    # Run minimap2 mapping
    echo "  Mapping reads with minimap2..."
    minimap2 -ax sr -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
        $REF_GENOME $R1 $R2 | \
        samtools view -@ $THREADS -bS - | \
        samtools sort -@ $THREADS -o $OUT_BAM

    # Index BAM file
    echo "  Indexing BAM file..."
    samtools index -@ $THREADS $OUT_BAM

    # Mark duplicates
    echo "  Marking and removing duplicates..."
    samtools markdup -r -@ $THREADS $OUT_BAM $OUTPUT_DIR/${sample}.dedup.bam
    samtools index -@ $THREADS $OUTPUT_DIR/${sample}.dedup.bam
    
    # Generate mapping statistics
    echo "  Generating mapping statistics..."
    samtools flagstat -@ $THREADS $OUTPUT_DIR/${sample}.dedup.bam > $OUTPUT_DIR/flagstats/${sample}.flagstat.txt
    
    # Generate coverage statistics
    echo "  Calculating coverage statistics..."
    samtools coverage -@ $THREADS $OUTPUT_DIR/${sample}.dedup.bam > $OUTPUT_DIR/coverage_stats/${sample}.coverage.txt
    
    # Calculate mean coverage
    mean_cov=$(samtools depth $OUTPUT_DIR/${sample}.dedup.bam | awk '{sum+=$3} END {print sum/NR}')
    echo "  Mean coverage depth: $mean_cov"
    
    echo "  Completed processing of $sample"
    echo "----------------------------------------"
}

# Process all samples in parallel
echo "Detecting samples..."
samples=$(ls $READS_DIR/*_R1_minimal.fq.gz | sed 's/.*\///' | sed 's/_R1_minimal.fq.gz//')
echo "Found samples: $samples"
echo "----------------------------------------"

# Process 1 sample at a time since we're using more threads per sample
for sample in $samples; do
    # Process this sample in the foreground since we're using 16 threads
    process_sample $sample
done

# Create a summary report
echo "Creating summary report..."
echo "Minimap2 Mapping Summary" > $OUTPUT_DIR/mapping_summary.txt
echo "======================" >> $OUTPUT_DIR/mapping_summary.txt
echo "Date: $(date)" >> $OUTPUT_DIR/mapping_summary.txt
echo "Reference genome: $REF_GENOME" >> $OUTPUT_DIR/mapping_summary.txt
echo "Number of samples processed: $(echo $samples | wc -w)" >> $OUTPUT_DIR/mapping_summary.txt
echo "Threads used per sample: $THREADS" >> $OUTPUT_DIR/mapping_summary.txt
echo "" >> $OUTPUT_DIR/mapping_summary.txt
echo "Mapping Statistics:" >> $OUTPUT_DIR/mapping_summary.txt
echo "----------------" >> $OUTPUT_DIR/mapping_summary.txt

for sample in $samples; do
    echo "Sample: $sample" >> $OUTPUT_DIR/mapping_summary.txt
    
    # Extract key metrics
    total=$(grep "in total" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | cut -d' ' -f1)
    mapped=$(grep "mapped (" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | head -1 | cut -d' ' -f1)
    mapping_rate=$(grep "mapped (" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | head -1 | grep -o "[0-9.]*%" | head -1)
    properly_paired=$(grep "properly paired" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | cut -d' ' -f1)
    proper_rate=$(grep "properly paired" $OUTPUT_DIR/flagstats/${sample}.flagstat.txt | grep -o "[0-9.]*%" | head -1)
    
    echo "  Total reads: $total" >> $OUTPUT_DIR/mapping_summary.txt
    echo "  Mapped reads: $mapped ($mapping_rate)" >> $OUTPUT_DIR/mapping_summary.txt
    echo "  Properly paired: $properly_paired ($proper_rate)" >> $OUTPUT_DIR/mapping_summary.txt
    echo "" >> $OUTPUT_DIR/mapping_summary.txt
done

# Calculate and record runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))

echo "Runtime: ${hours}h ${minutes}m ${seconds}s" >> $OUTPUT_DIR/mapping_summary.txt
echo "" >> $OUTPUT_DIR/mapping_summary.txt
echo "BAM files are in: $OUTPUT_DIR" >> $OUTPUT_DIR/mapping_summary.txt
echo "Flagstat files are in: $OUTPUT_DIR/flagstats" >> $OUTPUT_DIR/mapping_summary.txt
echo "Coverage statistics are in: $OUTPUT_DIR/coverage_stats" >> $OUTPUT_DIR/mapping_summary.txt

echo "Mapping process complete!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "Summary report is available at: $OUTPUT_DIR/mapping_summary.txt"
