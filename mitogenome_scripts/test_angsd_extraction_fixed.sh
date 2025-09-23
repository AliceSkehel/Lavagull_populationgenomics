#!/bin/bash

# ANGSD-based consensus - targeted to mtDNA only
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"

mkdir -p angsd_test_results

echo "Testing ANGSD method on mtDNA region only..."

# Test samples: LVGU_10, LVGU_1, LVGU_60, LVGU_39, LVGU_37
for sample in LVGU_10 LVGU_1 LVGU_60 LVGU_39 LVGU_37; do
    bam_file="illumina_bams/${sample}_Illumina.bam"
    
    if [ -f "$bam_file" ]; then
        echo "Processing $sample with ANGSD..."
        
        # Run ANGSD with region specification
        angsd -i $bam_file \
              -r OZ118781.1:1-16792 \
              -doFasta 2 \
              -doCounts 1 \
              -setMinDepth 4 \
              -minMapQ 30 \
              -out angsd_test_results/${sample}_angsd
        
        # Check if output was created and fix header
        if [ -f "angsd_test_results/${sample}_angsd.fa" ]; then
            sed -i "1s/.*/>$sample/" angsd_test_results/${sample}_angsd.fa
            echo "  Completed: $sample"
        else
            echo "  Warning: No FASTA output for $sample"
            ls angsd_test_results/${sample}_angsd*
        fi
    fi
done
