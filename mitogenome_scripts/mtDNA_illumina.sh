#!/bin/bash

# Simple Illumina consensus generation
mkdir -p ~/Sequencing_Combined/mtdna_only/consensus

for bam in ~/Sequencing_Combined/mtdna_only/LVGU_*_mtdna.bam; do
   sample=$(basename $bam _mtdna.bam)
   samtools consensus $bam > ~/Sequencing_Combined/mtdna_only/consensus/${sample}_consensus.fasta
   sed -i "1s/.*/>$sample/" ~/Sequencing_Combined/mtdna_only/consensus/${sample}_consensus.fasta
done

cat ~/Sequencing_Combined/mtdna_only/consensus/*_consensus.fasta > ~/Sequencing_Combined/mtdna_only/consensus/all_consensus.fasta

mafft --auto ~/Sequencing_Combined/mtdna_only/consensus/all_consensus.fasta > ~/Sequencing_Combined/mtdna_only/consensus/aligned.fasta

echo "Done!"
