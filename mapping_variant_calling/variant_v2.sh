#!/bin/bash

# Test ANGSD variant calling with 1 sample
BAM_DIR=~/Sequencing_Round2/complete_mapping
OUTPUT_DIR=~/Sequencing_Round2/variant_calling_test
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna

mkdir -p $OUTPUT_DIR

echo "Creating BAM list for test sample..."
echo "$BAM_DIR/LVGU_1.bam" > $OUTPUT_DIR/bam_list.txt

echo "Testing ANGSD variant calling on LVGU_1..."
angsd \
    -b $OUTPUT_DIR/bam_list.txt \
    -ref $REFERENCE \
    -out $OUTPUT_DIR/test_variants \
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
if [ -f "$OUTPUT_DIR/test_variants.bcf" ]; then
    echo "Converting to VCF..."
    bcftools view $OUTPUT_DIR/test_variants.bcf -Oz -o $OUTPUT_DIR/test_variants.vcf.gz
    bcftools index $OUTPUT_DIR/test_variants.vcf.gz
    echo "Success! Test VCF file: $OUTPUT_DIR/test_variants.vcf.gz"
    
    # Show some stats
    echo "Number of variants found: $(bcftools view -H $OUTPUT_DIR/test_variants.vcf.gz | wc -l)"
else
    echo "Error: ANGSD failed to create BCF file"
    echo "Check the ANGSD output above for errors"
fi

