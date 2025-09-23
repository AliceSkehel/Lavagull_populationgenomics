#!/bin/bash

# Test reference-based extraction method for Illumina data
# This method uses the reference sequence as template and only incorporates high-confidence variants
# Should eliminate soft-clipping artifacts causing alignment issues

REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"
REGION="OZ118781.1:1-16792"

mkdir -p test_reference_extraction

echo "Starting reference-based extraction test for Illumina samples..."
echo "Reference: $REFERENCE"
echo "Region: $REGION"

# Extract the reference template for this region
echo "Extracting reference template..."
samtools faidx $REFERENCE $REGION > test_reference_extraction/mtdna_template.fasta

# Test with a few Illumina samples
for sample in LVGU_1 LVGU_10 LVGU_20; do
    if [ -f "illumina_bams/${sample}_Illumina.bam" ]; then
        echo "Processing $sample with reference-based method..."
        
        # Generate VCF of high-confidence variants only
        bcftools mpileup -f $REFERENCE -r OZ118781.1:1-16792 \
            illumina_bams/${sample}_Illumina.bam | \
        bcftools call -mv --ploidy 1 | \
        bcftools filter -i 'QUAL>=30 && DP>=10' > test_reference_extraction/${sample}_variants.vcf
        
        # Apply only high-confidence variants to reference template
        bcftools consensus -f $REFERENCE test_reference_extraction/${sample}_variants.vcf \
            -s - > test_reference_extraction/${sample}_reference_based.fasta
        
        # Fix header
        sed -i "1s/.*/>$sample/" test_reference_extraction/${sample}_reference_based.fasta
        
        echo "  Completed: $sample"
    else
        echo "  Warning: BAM file not found for $sample"
    fi
done

echo "Test extraction complete!"
echo "Output files in: test_reference_extraction/"
echo "Compare these sequences to see if alignment issues are resolved."
