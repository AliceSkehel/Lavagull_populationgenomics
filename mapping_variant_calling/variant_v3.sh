#!/bin/bash

# ANGSD variant calling for all samples
BAM_DIR=~/Sequencing_Round2/complete_mapping
OUTPUT_DIR=~/Sequencing_Round2/variant_calling
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna

mkdir -p $OUTPUT_DIR

echo "Creating BAM list for all samples..."
ls $BAM_DIR/LVGU_*.bam > $OUTPUT_DIR/bam_list.txt

echo "Running ANGSD variant calling on all $(wc -l < $OUTPUT_DIR/bam_list.txt) samples..."
angsd \
    -b $OUTPUT_DIR/bam_list.txt \
    -ref $REFERENCE \
    -out $OUTPUT_DIR/variants \
    -minMapQ 30 \
    -minQ 20 \
    -minMaf 0.05 \
    -GL 1 \
    -doMaf 1 \
    -doMajorMinor 1 \
    -doPost 1 \
    -doGeno 1 \
    -doBcf 1 \
    -SNP_pval 1e-6 \
    -P 8

# Check if BCF was created
if [ -f "$OUTPUT_DIR/variants.bcf" ]; then
    echo "Converting to VCF..."
    bcftools view $OUTPUT_DIR/variants.bcf -Oz -o $OUTPUT_DIR/variants.vcf.gz
    bcftools index $OUTPUT_DIR/variants.vcf.gz
    echo "Success! VCF file: $OUTPUT_DIR/variants.vcf.gz"
    
    # Show some stats
    echo "Number of variants found: $(bcftools view -H $OUTPUT_DIR/variants.vcf.gz | wc -l)"
    echo "Number of samples: $(bcftools query -l $OUTPUT_DIR/variants.vcf.gz | wc -l)"
else
    echo "Error: ANGSD failed to create BCF file"
    echo "Check the ANGSD output above for errors"
fi
