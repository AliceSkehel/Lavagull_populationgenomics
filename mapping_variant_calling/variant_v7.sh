#!/bin/bash

# Fast ANGSD - Illumina Only (32 samples)
WORK_DIR=~/Sequencing_Combined/angsd_variant_calling
mkdir -p $WORK_DIR

# ANGSD with illumina-only samples
angsd -bam ~/illumina_only_bam_list.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_illumina_only \
      -GL 1 \
      -doGlf 2 \
      -doMajorMinor 1 \
      -doMaf 2 \
      -doGeno 3 \
      -doPost 1 \
      -uniqueOnly 1 \
      -remove_bads 1 \
      -only_proper_pairs 1 \
      -trim 0 \
      -baq 1 \
      -minInd 25 \
      -doCounts 1 \
      -setMinDepth 120 \
      -setMaxDepth 2534 \
      -setMinDepthInd 3 \
      -setMaxDepthInd 30 \
      -rf ~/regions_no_mtdna.txt \
      -nThreads 12

echo "Done!"
