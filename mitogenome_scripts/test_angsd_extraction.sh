#!/bin/bash

# ANGSD-based consensus for specific test samples
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"

mkdir -p angsd_test_results

echo "Testing ANGSD method on selected samples..."

# Test samples: LVGU_10, LVGU_1, LVGU_60, LVGU_39, LVGU_37
for sample in LVGU_10 LVGU_1 LVGU_60 LVGU_39 LVGU_37; do
    bam_file="illumina_bams/${sample}_Illumina.bam"
    
    if [ -f "$bam_file" ]; then
        echo "Processing $sample with ANGSD..."
        
        # Extract mtDNA region only
        samtools view -b $bam_file OZ118781.1:1-16792 > angsd_test_results/${sample}_mtdna_only.bam
        samtools index angsd_test_results/${sample}_mtdna_only.bam
        
        # Run ANGSD consensus
        angsd -i angsd_test_results/${sample}_mtdna_only.bam \
              -doFasta 2 \
              -doCounts 1 \
              -setMinDepth 4 \
              -minMapQ 30 \
              -explode 1 \
              -out angsd_test_results/${sample}_angsd
        
        # Fix header if FASTA was created
        if [ -f "angsd_test_results/${sample}_angsd.fa" ]; then
            sed -i "1s/.*/>$sample/" angsd_test_results/${sample}_angsd.fa
            echo "  Completed: $sample"
        else
            echo "  Warning: No FASTA output for $sample"
        fi
    else
        echo "  Warning: BAM file not found for $sample"
    fi
done

echo "ANGSD test complete!"
echo "Results in: angsd_test_results/"
