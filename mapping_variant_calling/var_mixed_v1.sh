#!/bin/bash

# ANGSD Variant Calling - Lava Gulls (VCF output)
WORK_DIR=~/Sequencing_Combined/angsd_variant_calling
mkdir -p $WORK_DIR

# Create BAM list
ls ~/Sequencing_Combined/complete_mapping/LVGU_*.bam > $WORK_DIR/all_bam_list.txt

echo "Starting ANGSD..."

# Run ANGSD for VCF output
angsd -bam $WORK_DIR/all_bam_list.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_all33 \
      -GL 1 \
      -doBcf 1 \
      -doMajorMinor 1 \
      -doMaf 2 \
      -doGeno 3 \
      -doPost 1 \
      -postCutoff 0.95 \
      -uniqueOnly 1 \
      -remove_bads 1 \
      -only_proper_pairs 1 \
      -trim 0 \
      -baq 1 \
      -minMapQ 20 \
      -minQ 20 \
      -minInd 25 \
      -minMaf 0.05 \
      -SNP_pval 1e-6 \
      -doCounts 1 \
      -setMinDepth 75 \
      -setMaxDepth 4356 \
      -skipTriallelic 1 \
      -rf ~/regions_no_mtdna.txt \
      -nThreads 8

echo "Done!"

