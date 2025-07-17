#!/bin/bash

# Simple ONT mtDNA Consensus Generation
mkdir -p ~/mitogenome/consensus

echo "Generating ONT mtDNA consensus sequences..."

# Good quality samples
GOOD_SAMPLES=(
    "LVGU_001" "LVGU_004" "LVGU_005" "LVGU_006" "LVGU_010" 
    "LVGU_012" "LVGU_014" "LVGU_018" "LVGU_020" "LVGU_023" 
    "LVGU_030" "LVGU_034" "LVGU_039" "LVGU_041" "LVGU_042"
)

# Generate consensus for each sample
for sample in "${GOOD_SAMPLES[@]}"; do
    echo "Processing $sample..."
    
    samtools consensus ~/mitogenome/mtdna_only/${sample}_mtDNA_ONT.bam > ~/mitogenome/consensus/${sample}_consensus.fasta
    sed -i "1s/.*/>$sample/" ~/mitogenome/consensus/${sample}_consensus.fasta
    
    echo "  ✅ Done: ${sample}_consensus.fasta"
done

# Combine all consensus sequences
cat ~/mitogenome/consensus/*_consensus.fasta > ~/mitogenome/consensus/all_consensus.fasta

echo ""
echo "🎉 Consensus generation complete!"
echo "Combined file: ~/mitogenome/consensus/all_consensus.fasta"
echo ""

# Align sequences
echo "Aligning consensus sequences..."
mafft --auto ~/mitogenome/consensus/all_consensus.fasta > ~/mitogenome/consensus/aligned_mtDNA.fasta

echo "✅ Aligned sequences: ~/mitogenome/consensus/aligned_mtDNA.fasta"
echo "Ready for phylogenetic analysis!"
