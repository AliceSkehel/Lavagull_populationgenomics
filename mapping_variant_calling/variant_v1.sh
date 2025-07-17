#!/bin/bash

# Simple ANGSD variant calling
BAM_DIR=~/Sequencing_Round2/complete_mapping
OUTPUT_DIR=~/Sequencing_Round2/variant_calling
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna

mkdir -p $OUTPUT_DIR

echo "Creating BAM list..."
ls $BAM_DIR/LVGU_*.bam > $OUTPUT_DIR/bam_list.txt

echo "Running ANGSD variant calling..."
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
    -doVcf 1 \
    -doBcf 1 \
    -SNP_pval 1e-6 \
    -P 8

echo "Converting to VCF..."
bcftools view $OUTPUT_DIR/variants.bcf -Oz -o $OUTPUT_DIR/variants.vcf.gz
bcftools index $OUTPUT_DIR/variants.vcf.gz

echo "Done! VCF file: $OUTPUT_DIR/variants.vcf.gz"
