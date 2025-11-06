#!/bin/bash

GENOME_SIZE=1.31

echo -e "Sample\tSequencing_Round\tRaw_R1_Gbp\tRaw_R2_Gbp\tRaw_Total_Gbp\tRaw_Coverage"

# Process Round 2 samples
echo "Processing Round 2 samples..." >&2
for file in /home/askehel/Sequencing_Round2/fastq_Illumina_June/01.RawData/*_1.fq.gz; do
    if [ -f "$file" ]; then
        sample=$(basename $file _1.fq.gz)
        echo "Processing $sample..." >&2
        
        raw_r1=$(zcat ${file} 2>/dev/null | awk 'NR%4==2 {sum+=length($0)} END {print sum/1e9}')
        raw_r2=$(zcat ${file/_1.fq.gz/_2.fq.gz} 2>/dev/null | awk 'NR%4==2 {sum+=length($0)} END {print sum/1e9}')
        raw_total=$(echo "$raw_r1 + $raw_r2" | bc -l)
        raw_cov=$(echo "scale=2; $raw_total / $GENOME_SIZE" | bc -l)
        
        echo -e "$sample\tRound2\t$raw_r1\t$raw_r2\t$raw_total\t$raw_cov"
    fi
done

# Process Round 1 samples
echo "Processing Round 1 samples..." >&2
for file in /home/askehel/Sequencing_Round1/*_L1_1.fq.gz; do
    if [ -f "$file" ]; then
        sample=$(basename $file _L1_1.fq.gz)
        echo "Processing $sample..." >&2
        
        raw_r1=$(zcat ${file} 2>/dev/null | awk 'NR%4==2 {sum+=length($0)} END {print sum/1e9}')
        raw_r2=$(zcat ${file/_L1_1.fq.gz/_L1_2.fq.gz} 2>/dev/null | awk 'NR%4==2 {sum+=length($0)} END {print sum/1e9}')
        raw_total=$(echo "$raw_r1 + $raw_r2" | bc -l)
        raw_cov=$(echo "scale=2; $raw_total / $GENOME_SIZE" | bc -l)
        
        echo -e "$sample\tRound1\t$raw_r1\t$raw_r2\t$raw_total\t$raw_cov"
    fi
done
