#!/bin/bash

# ANGSD processing for ancient DNA samples
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"

mkdir -p ancient_angsd_results

echo "Processing ancient DNA samples with ANGSD..."

# Process both ancient samples with more permissive parameters
for sample in LVGU_473 LVGU_476; do
    bam_file="Sequencing_Round2/ancient_mapping_ont/${sample}.YellowLeggedGull.bam"
    
    if [ -f "$bam_file" ]; then
        echo "Processing ancient sample: $sample"
        
        # More permissive parameters for ancient DNA
        angsd -i $bam_file \
              -r OZ118781.1:1-16792 \
              -doFasta 2 \
              -doCounts 1 \
              -setMinDepth 2 \
              -minMapQ 20 \
              -out ancient_angsd_results/${sample}_ancient_angsd
        
        # Decompress and fix header
        if [ -f "ancient_angsd_results/${sample}_ancient_angsd.fa.gz" ]; then
            gunzip ancient_angsd_results/${sample}_ancient_angsd.fa.gz
            sed -i "1s/.*/>$sample/" ancient_angsd_results/${sample}_ancient_angsd.fa
            echo "  Completed: $sample"
        fi
    fi
done

# Combine ancient sequences
cat ancient_angsd_results/*.fa > ancient_all_angsd.fasta
echo "Ancient DNA processing complete!"
