#!/bin/bash

REFGEN=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna
OUTDIR=~/Sequencing_Combined/angsd_variant_calling_all35/heterozygosity_analysis
BAMLIST=~/all_samples_bam_list_no_LVGU60.txt

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

while read BAM; do
    SAMPLE=$(basename $BAM .bam | sed 's/.YellowLeggedGull_ONT//g')
    echo "Processing $SAMPLE..."
    
    # Calculate SAF (removed -setMinDepthInd and -setMaxDepthInd)
    angsd -i $BAM -ref $REFGEN -anc $REFGEN \
          -doSaf 1 -GL 1 \
          -minMapQ 30 -minQ 20 \
          -C 50 \
          -out $OUTDIR/${SAMPLE}_saf
    
    # Calculate SFS
    realSFS $OUTDIR/${SAMPLE}_saf.saf.idx > $OUTDIR/${SAMPLE}.sfs
    
done < $BAMLIST

echo "All heterozygosity calculations complete!"
