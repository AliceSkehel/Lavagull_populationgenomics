#!/bin/bash

# Reference-based extraction for ALL Illumina samples
# Uses mtDNA template to avoid whole-genome processing

REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"
REGION="OZ118781.1:1-16792"

mkdir -p illumina_reference_based

echo "Starting reference-based extraction for all Illumina samples..."

# Extract mtDNA reference template
echo "Creating mtDNA template..."
samtools faidx $REFERENCE $REGION > illumina_reference_based/mtdna_template.fasta

# Process all Illumina samples
for bam in illumina_bams/*_Illumina.bam; do
    if [ -f "$bam" ]; then
        sample=$(basename $bam _Illumina.bam)
        echo "Processing $sample..."
        
        # Generate high-confidence variants for mtDNA region only
        bcftools mpileup -f $REFERENCE -r $REGION $bam | \
        bcftools call -mv --ploidy 1 | \
        bcftools filter -i 'QUAL>=30 && DP>=10' | \
        bgzip > illumina_reference_based/${sample}_variants.vcf.gz
        
        tabix -p vcf illumina_reference_based/${sample}_variants.vcf.gz
        
        # Apply variants to mtDNA template only (not whole genome)
        bcftools consensus -f illumina_reference_based/mtdna_template.fasta \
            illumina_reference_based/${sample}_variants.vcf.gz > \
            illumina_reference_based/${sample}_reference_based.fasta
        
        # Fix header
        sed -i "1s/.*/>$sample/" illumina_reference_based/${sample}_reference_based.fasta
        
        echo "  Completed: $sample"
    fi
done

echo "All Illumina samples processed!"
echo "Files in: illumina_reference_based/"
