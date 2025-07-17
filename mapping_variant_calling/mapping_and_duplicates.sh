#!/bin/bash
# map_lagu_samples_fast.sh - Optimized script to map LaGu samples to European Herring Gull reference

# Reference genome path
REF="/home/askehel/Doves_run/genome_data/GCA_964417175.1_bLarArg3.hap1.1_genomic.fna"

# Output directory
OUT_DIR="/home/askehel/Doves_run/mapping_results"

# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# Number of threads per sample
THREADS=8  # Adjust based on your system resources
MEMORY="4G"  # Memory per thread for sorting

# Function to process a sample
process_sample() {
    sample_dir=$1
    sample_name=$(basename $sample_dir)
    
    # Create a log file
    log_file="$OUT_DIR/${sample_name}.log"
    echo "Processing sample: $sample_name (Logging to $log_file)" | tee $log_file
    
    # Find fastq files
    read1_files=$(find /home/askehel/Doves_run/$sample_dir -name "*_1_fastp_clean.fq.gz" -o -name "*_1.fq.gz" | sort)
    read2_files=$(find /home/askehel/Doves_run/$sample_dir -name "*_2_fastp_clean.fq.gz" -o -name "*_2.fq.gz" | sort)
    
    # If no files found, try alternative patterns
    if [ -z "$read1_files" ]; then
        read1_files=$(find /home/askehel/Doves_run/$sample_dir -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" | sort)
        read2_files=$(find /home/askehel/Doves_run/$sample_dir -name "*_R2_*.fastq.gz" -o -name "*_R2_*.fq.gz" | sort)
    fi
    
    # Check if fastq files were found
    if [ -z "$read1_files" ] || [ -z "$read2_files" ]; then
        echo "Could not find FASTQ files in /home/askehel/Doves_run/$sample_dir. Skipping." | tee -a $log_file
        return
    fi
    
    # Convert to lists for BWA
    read1_list=$(echo "$read1_files" | tr '\n' ' ')
    read2_list=$(echo "$read2_files" | tr '\n' ' ')
    
    echo "Starting mapping with $THREADS threads..." | tee -a $log_file
    
    # Streamlined pipeline with pipes
    echo "Running streamlined mapping pipeline..." | tee -a $log_file
    bwa mem -t $THREADS $REF $read1_list $read2_list | \
        samtools sort -@ $THREADS -m $MEMORY -o $OUT_DIR/${sample_name}.sorted.bam 2>> $log_file
    
    # Index the BAM file
    echo "Indexing BAM file..." | tee -a $log_file
    samtools index $OUT_DIR/${sample_name}.sorted.bam 2>> $log_file
    
    # Optimized duplicate marking process
    echo "Processing duplicates..." | tee -a $log_file
    samtools sort -n -@ $THREADS -m $MEMORY -o $OUT_DIR/${sample_name}.namesorted.bam $OUT_DIR/${sample_name}.sorted.bam 2>> $log_file
    samtools fixmate -m -@ $THREADS $OUT_DIR/${sample_name}.namesorted.bam $OUT_DIR/${sample_name}.fixmate.bam 2>> $log_file
    samtools sort -@ $THREADS -m $MEMORY -o $OUT_DIR/${sample_name}.positionsorted.bam $OUT_DIR/${sample_name}.fixmate.bam 2>> $log_file
    samtools markdup -@ $THREADS $OUT_DIR/${sample_name}.positionsorted.bam $OUT_DIR/${sample_name}.marked.bam 2>> $log_file
    samtools index $OUT_DIR/${sample_name}.marked.bam 2>> $log_file
    
    # Get duplicate statistics
    echo "Getting duplicate statistics..." | tee -a $log_file
    samtools flagstat -@ $THREADS $OUT_DIR/${sample_name}.marked.bam | grep -A 2 "duplicates" | tee -a $log_file
    
    # Remove duplicates
    echo "Removing duplicates..." | tee -a $log_file
    samtools markdup -r -@ $THREADS $OUT_DIR/${sample_name}.positionsorted.bam $OUT_DIR/${sample_name}.dedup.bam 2>> $log_file
    samtools index $OUT_DIR/${sample_name}.dedup.bam 2>> $log_file
    
    # Parallel statistics calculation
    echo "Calculating statistics..." | tee -a $log_file
    samtools flagstat -@ $THREADS $OUT_DIR/${sample_name}.dedup.bam > $OUT_DIR/${sample_name}.flagstat.txt &
    samtools coverage -@ $THREADS $OUT_DIR/${sample_name}.dedup.bam > $OUT_DIR/${sample_name}.coverage.txt &
    
    # Wait for parallel stats to complete
    wait
    
    # Calculate average depth
    echo "Calculating depth (this might take a while)..." | tee -a $log_file
    samtools depth -a $OUT_DIR/${sample_name}.dedup.bam | \
        awk '{sum+=$3} END {print "Average depth = " sum/NR "X"}' > $OUT_DIR/${sample_name}.depth.txt
    cat $OUT_DIR/${sample_name}.depth.txt | tee -a $log_file
    
    # Clean up intermediate files
    echo "Cleaning up intermediate files..." | tee -a $log_file
    rm $OUT_DIR/${sample_name}.namesorted.bam
    rm $OUT_DIR/${sample_name}.fixmate.bam
    
    echo "Completed processing for $sample_name" | tee -a $log_file
    echo "============================================================" | tee -a $log_file
}

# Process samples in parallel
run_parallel() {
    local max_concurrent=2  # Number of samples to process simultaneously
    
    echo "Starting mapping of LaGu samples in parallel (max $max_concurrent at once)"
    echo "Using $THREADS threads per sample"
    
    for sample in LaGu_Workshop05 LaGu_Workshop06 LaGu_Workshop08 LaGu_Workshop09A LaGu_Workshop09B; do
        # Wait if we already have max jobs running
        while [ $(jobs -r | wc -l) -ge $max_concurrent ]; do
            sleep 10
        done
        
        echo "Launching process for $sample"
        process_sample "$sample" &
        
        # Brief pause to prevent resource contention on startup
        sleep 5
    done
    
    # Wait for all background jobs to complete
    wait
    echo "All sample processing complete!"
}

# Run the parallel processing function
run_parallel

# Create a summary file after all processing is complete
echo "Creating summary file..."
echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapping_Rate\tDuplicate_Rate\tAverage_Depth" > $OUT_DIR/lagu_mapping_summary.txt

for sample in LaGu_Workshop04 LaGu_Workshop05 LaGu_Workshop06 LaGu_Workshop08 LaGu_Workshop09A LaGu_Workshop09B; do
    # Skip if files don't exist
    if [ ! -f "$OUT_DIR/${sample}.flagstat.txt" ]; then
        echo "$sample: Skipping - files not found"
        continue
    fi
    
    # Extract statistics
    total=$(grep "in total" $OUT_DIR/${sample}.flagstat.txt | awk '{print $1}')
    mapped=$(grep "mapped (" $OUT_DIR/${sample}.flagstat.txt | awk '{print $1}')
    map_rate=$(grep "mapped (" $OUT_DIR/${sample}.flagstat.txt | awk -F "(" '{print $2}' | awk -F ":" '{print $1}')
    
    # Get duplicate rate from marked BAM
    if [ -f "$OUT_DIR/${sample}.marked.bam" ]; then
        dup_count=$(samtools flagstat $OUT_DIR/${sample}.marked.bam | grep "duplicates" | head -1 | awk '{print $1}')
        dup_rate=$(echo "scale=2; 100*$dup_count/$total" | bc)
    else
        dup_rate="NA"
    fi
    
    # Get average depth
    if [ -f "$OUT_DIR/${sample}.depth.txt" ]; then
        avg_depth=$(cat $OUT_DIR/${sample}.depth.txt | awk '{print $4}' | sed 's/X//')
    else
        avg_depth="NA"
    fi
    
    echo -e "$sample\t$total\t$mapped\t$map_rate\t$dup_rate\t$avg_depth" >> $OUT_DIR/lagu_mapping_summary.txt
done

echo "Summary created at $OUT_DIR/lagu_mapping_summary.txt"
echo "All LaGu samples processed successfully!"
