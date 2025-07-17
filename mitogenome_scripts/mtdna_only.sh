#!/bin/bash

# Create directory for mtDNA-only files
mkdir -p ~/mitogenome/mtdna_only

# Extract mtDNA reads from all samples
for bam in ~/mitogenome/mapped_bams/*_ONT.bam; do
    sample=$(basename $bam _ONT.bam)
    echo "Extracting mtDNA for $sample..."
    
    samtools view -b $bam "OZ118781.1" > ~/mitogenome/mtdna_only/${sample}_mtDNA_ONT.bam
    samtools index ~/mitogenome/mtdna_only/${sample}_mtDNA_ONT.bam
    
    # Quick read count
    reads=$(samtools view -c ~/mitogenome/mtdna_only/${sample}_mtDNA_ONT.bam)
    echo "  $sample: $reads mtDNA reads"
done
