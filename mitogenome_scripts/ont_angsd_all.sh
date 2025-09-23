#!/bin/bash

# ANGSD-based consensus for ALL corrected ONT samples
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"

mkdir -p ont_angsd_all

echo "Processing all corrected ONT samples with ANGSD..."

# Process all corrected ONT samples
for bam in /home/askehel/mitogenome/corrected_mapping/*_ONT_corrected.bam; do
    if [ -f "$bam" ]; then
        sample=$(basename $bam _ONT_corrected.bam)
        echo "Processing $sample with ANGSD..."
        
        # Run ANGSD with same parameters as Illumina
        angsd -i $bam \
              -r OZ118781.1:1-16792 \
              -doFasta 2 \
              -doCounts 1 \
              -setMinDepth 4 \
              -minMapQ 30 \
              -out ont_angsd_all/${sample}_angsd
        
        # Decompress if compressed
        if [ -f "ont_angsd_all/${sample}_angsd.fa.gz" ]; then
            gunzip ont_angsd_all/${sample}_angsd.fa.gz
        fi
        
        # Fix header
        if [ -f "ont_angsd_all/${sample}_angsd.fa" ]; then
            sed -i "1s/.*/>$sample/" ont_angsd_all/${sample}_angsd.fa
            echo "  Completed: $sample"
        else
            echo "  Warning: No FASTA output for $sample"
        fi
    fi
done

# Combine all ONT sequences
cat ont_angsd_all/*.fa > ont_all_angsd.fasta

echo "All ONT samples processed with ANGSD!"
echo "Combined file: ont_all_angsd.fasta"
echo "Individual files in: ont_angsd_all/"
