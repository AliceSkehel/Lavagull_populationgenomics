#!/bin/bash

# ANGSD genotype likelihood analysis for subset samples
BAM_DIR=~/Sequencing_Round2/complete_mapping/subset
OUTPUT_DIR=~/Sequencing_Round2/variant_calling
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna

mkdir -p $OUTPUT_DIR

echo "Creating BAM list for subset samples..."
ls $BAM_DIR/LVGU_*.bam > $OUTPUT_DIR/bam_list.txt

echo "Running ANGSD genotype likelihood analysis on $(wc -l < $OUTPUT_DIR/bam_list.txt) samples..."
angsd \
    -b $OUTPUT_DIR/bam_list.txt \
    -ref $REFERENCE \
    -out $OUTPUT_DIR/variants \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -trim 0 \
    -baq 1 \
    -minMapQ 30 \
    -minQ 20 \
    -minInd 16 \
    -setMinDepth 115 \
    -setMaxDepth 845 \
    -minMaf 0.05 \
    -GL 1 \
    -doMaf 1 \
    -doMajorMinor 1 \
    -doPost 1 \
    -doGlf 2 \
    -doDepth 1 \
    -doCounts 1 \
    -SNP_pval 1e-6 \
    -P 8

# Check if genotype likelihood files were created
if [ -f "$OUTPUT_DIR/variants.beagle.gz" ]; then
    echo "Success! Genotype likelihood file: $OUTPUT_DIR/variants.beagle.gz"
    echo "MAF file: $OUTPUT_DIR/variants.mafs.gz"
    
    # Show some stats
    echo "Number of sites: $(zcat $OUTPUT_DIR/variants.mafs.gz | tail -n +2 | wc -l)"
    echo "Number of samples: $(zcat $OUTPUT_DIR/variants.beagle.gz | head -1 | tr '\t' '\n' | tail -n +4 | wc -l | awk '{print $1/3}')"
else
    echo "Error: ANGSD failed to create genotype likelihood files"
    echo "Check the ANGSD output above for errors"
fi
