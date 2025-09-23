#!/bin/bash

# ANGSD-based consensus for ALL Illumina samples
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"

mkdir -p illumina_angsd_all

echo "Processing all Illumina samples with ANGSD..."

# Process all Illumina samples
for bam in illumina_bams/*_Illumina.bam; do
    if [ -f "$bam" ]; then
        sample=$(basename $bam _Illumina.bam)
        echo "Processing $sample with ANGSD..."
        
        # Run ANGSD with region specification
        angsd -i $bam \
              -r OZ118781.1:1-16792 \
              -doFasta 2 \
              -doCounts 1 \
              -setMinDepth 4 \
              -minMapQ 30 \
              -out illumina_angsd_all/${sample}_angsd
        
        # Decompress if compressed
        if [ -f "illumina_angsd_all/${sample}_angsd.fa.gz" ]; then
            gunzip illumina_angsd_all/${sample}_angsd.fa.gz
        fi
        
        # Fix header
        if [ -f "illumina_angsd_all/${sample}_angsd.fa" ]; then
            sed -i "1s/.*/>$sample/" illumina_angsd_all/${sample}_angsd.fa
            echo "  Completed: $sample"
        else
            echo "  Warning: No FASTA output for $sample"
        fi
    fi
done

# Combine all sequences
cat illumina_angsd_all/*.fa > illumina_all_angsd.fasta

echo "All Illumina samples processed with ANGSD!"
echo "Combined file: illumina_all_angsd.fasta"
echo "Individual files in: illumina_angsd_all/"
