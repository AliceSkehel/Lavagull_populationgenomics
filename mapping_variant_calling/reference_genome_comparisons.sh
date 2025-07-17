#!/bin/bash

# Define paths
CLEAN_READS=~/Sequencing_Round1/fastp_cleaned_reads
HERRING_GULL_REF=~/Doves_run/genome_data/GCA_964417175.1_bLarArg3.hap1.1_genomic.fna
YELLOW_GULL_REF=~/Doves_run/genome_data/yellow_legged_gull_genome.fna
OUTPUT_DIR=~/Sequencing_Round1/reference_comparison
COMPARISON_REPORT=$OUTPUT_DIR/reference_comparison_report.txt
THREADS=8 # Set number of threads consistently

# Create output directory
mkdir -p $OUTPUT_DIR

# Specify the two samples to compare
SAMPLES="LaGu_Workshop06 LaGu_Workshop08"

# Create header for comparison report
echo "Reference Genome Mapping Comparison" > $COMPARISON_REPORT
echo "=================================" >> $COMPARISON_REPORT
echo "" >> $COMPARISON_REPORT
date >> $COMPARISON_REPORT
echo "" >> $COMPARISON_REPORT
echo "Samples: LaGu_Workshop06 and LaGu_Workshop08" >> $COMPARISON_REPORT
echo "Threads used: $THREADS" >> $COMPARISON_REPORT
echo "" >> $COMPARISON_REPORT

# Function to process a sample with a specific reference
process_mapping() {
    local sample=$1
    local ref_type=$2
    local ref_path=$3
    local r1=$4
    local r2=$5
    
    echo "Mapping $sample to $ref_type reference..."
    
    # Map to reference
    bwa mem -t $THREADS -M -R "@RG\tID:${sample}_${ref_type}\tSM:${sample}\tPL:ILLUMINA" \
        $ref_path \
        $r1 $r2 \
        | samtools view -@ $THREADS -bS - \
        | samtools sort -@ $THREADS -o $OUTPUT_DIR/${sample}_${ref_type}.sorted.bam
    
    # Index BAM (use threads if available)
    samtools index -@ $THREADS $OUTPUT_DIR/${sample}_${ref_type}.sorted.bam
    
    # Calculate mapping stats (use threads if available)
    echo "Calculating statistics for $ref_type mapping..."
    samtools flagstat -@ $THREADS $OUTPUT_DIR/${sample}_${ref_type}.sorted.bam > $OUTPUT_DIR/${sample}_${ref_type}.flagstat.txt
    
    # Calculate mapping quality distribution
    samtools view -@ $THREADS -F 0x904 $OUTPUT_DIR/${sample}_${ref_type}.sorted.bam | \
        cut -f5 | sort -n | uniq -c > $OUTPUT_DIR/${sample}_${ref_type}.mapq.txt
    
    # Calculate coverage statistics
    samtools coverage -@ $THREADS $OUTPUT_DIR/${sample}_${ref_type}.sorted.bam -o $OUTPUT_DIR/${sample}_${ref_type}.coverage.txt
    
    echo "Completed $ref_type mapping for $sample"
}

# Process each sample
for sample in $SAMPLES; do
    echo "Processing sample: $sample"
    echo "Sample: $sample" >> $COMPARISON_REPORT
    echo "-----------------" >> $COMPARISON_REPORT
    
    # Define input files
    R1=$CLEAN_READS/${sample}_R1_cleaned.fq.gz
    R2=$CLEAN_READS/${sample}_R2_cleaned.fq.gz
    
    # Check if files exist
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Error: Input files for $sample not found. Skipping."
        echo "Error: Input files for $sample not found. Skipping." >> $COMPARISON_REPORT
        continue
    fi
    
    # Run mappings in parallel using background processes
    process_mapping "$sample" "herring" "$HERRING_GULL_REF" "$R1" "$R2" &
    process_mapping "$sample" "yellow" "$YELLOW_GULL_REF" "$R1" "$R2" &
    
    # Wait for both mapping processes to complete
    wait
    
    # Extract and compare key metrics
    HERRING_MAPPED=$(grep "mapped (" $OUTPUT_DIR/${sample}_herring.flagstat.txt | head -1 | cut -d' ' -f1)
    HERRING_TOTAL=$(grep "in total" $OUTPUT_DIR/${sample}_herring.flagstat.txt | cut -d' ' -f1)
    HERRING_PROPERLY_PAIRED=$(grep "properly paired" $OUTPUT_DIR/${sample}_herring.flagstat.txt | cut -d' ' -f1)
    
    YELLOW_MAPPED=$(grep "mapped (" $OUTPUT_DIR/${sample}_yellow.flagstat.txt | head -1 | cut -d' ' -f1)
    YELLOW_TOTAL=$(grep "in total" $OUTPUT_DIR/${sample}_yellow.flagstat.txt | cut -d' ' -f1)
    YELLOW_PROPERLY_PAIRED=$(grep "properly paired" $OUTPUT_DIR/${sample}_yellow.flagstat.txt | cut -d' ' -f1)
    
    # Calculate percentages
    HERRING_MAPPING_RATE=$(echo "scale=2; 100 * $HERRING_MAPPED / $HERRING_TOTAL" | bc)
    HERRING_PROPER_PAIR_RATE=$(echo "scale=2; 100 * $HERRING_PROPERLY_PAIRED / $HERRING_TOTAL" | bc)
    
    YELLOW_MAPPING_RATE=$(echo "scale=2; 100 * $YELLOW_MAPPED / $YELLOW_TOTAL" | bc)
    YELLOW_PROPER_PAIR_RATE=$(echo "scale=2; 100 * $YELLOW_PROPERLY_PAIRED / $YELLOW_TOTAL" | bc)
    
    # Write comparison to report
    echo "European Herring Gull Reference:" >> $COMPARISON_REPORT
    echo "  Total reads: $HERRING_TOTAL" >> $COMPARISON_REPORT
    echo "  Mapped reads: $HERRING_MAPPED (${HERRING_MAPPING_RATE}%)" >> $COMPARISON_REPORT
    echo "  Properly paired reads: $HERRING_PROPERLY_PAIRED (${HERRING_PROPER_PAIR_RATE}%)" >> $COMPARISON_REPORT
    echo "" >> $COMPARISON_REPORT
    
    echo "Yellow-legged Gull Reference:" >> $COMPARISON_REPORT
    echo "  Total reads: $YELLOW_TOTAL" >> $COMPARISON_REPORT
    echo "  Mapped reads: $YELLOW_MAPPED (${YELLOW_MAPPING_RATE}%)" >> $COMPARISON_REPORT
    echo "  Properly paired reads: $YELLOW_PROPERLY_PAIRED (${YELLOW_PROPER_PAIR_RATE}%)" >> $COMPARISON_REPORT
    echo "" >> $COMPARISON_REPORT
    
    # Simple recommendation
    if (( $(echo "$HERRING_MAPPING_RATE > $YELLOW_MAPPING_RATE" | bc -l) )); then
        BETTER="European Herring Gull"
        DIFF=$(echo "$HERRING_MAPPING_RATE - $YELLOW_MAPPING_RATE" | bc)
    else
        BETTER="Yellow-legged Gull"
        DIFF=$(echo "$YELLOW_MAPPING_RATE - $HERRING_MAPPING_RATE" | bc)
    fi
    
    echo "Recommendation for $sample: $BETTER reference performed better by ${DIFF}% mapping rate" >> $COMPARISON_REPORT
    echo "==========================================================" >> $COMPARISON_REPORT
    echo "" >> $COMPARISON_REPORT
    
    echo "Completed processing of $sample"
done

# Process both samples in parallel
echo "Processing both samples in parallel with 8 threads each..."

# Add overall recommendation based on both samples
echo "OVERALL RECOMMENDATION" >> $COMPARISON_REPORT
echo "======================" >> $COMPARISON_REPORT
echo "Check the individual sample results above to determine which reference" >> $COMPARISON_REPORT
echo "performed better consistently across both samples." >> $COMPARISON_REPORT
echo "" >> $COMPARISON_REPORT
echo "After determining the better reference, use it for mapping all samples." >> $COMPARISON_REPORT

echo "Comparison complete! Results available in $COMPARISON_REPORT"
