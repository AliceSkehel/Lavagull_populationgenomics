#!/bin/bash

# Samtools for all 35 samples
mkdir -p ~/qualimap_ont_filtered_samples

echo "Calculating coverage for all 35 samples......"
echo "Start time: $(date)"

for bam in $(cat ~/all_samples_bam_list.txt); do
    sample=$(basename $bam .bam)
    echo -n "$sample: "
    samtools depth $bam | awk '{sum+=$3; count++} END {print sum/count, "x"}'
done > ~/sample_coverage_summary.txt

echo "Done!"
echo "End time: $(date)"
cat ~/sample_coverage_summary.txt
