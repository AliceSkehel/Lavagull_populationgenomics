#!/bin/bash

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/modern_samples_only_bam_list.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

cd $WORK_DIR

echo "Calculating thetas for Group 3 (modern samples only, n=33)..."
echo "Start time: $(date)"

realSFS lava_gulls_saf_group3.saf.idx > lava_gulls_group3_modern.sfs

angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -rf ~/regions_group3.txt \
      -doThetas 1 -doSaf 1 \
      -pest lava_gulls_group3_modern.sfs \
      -GL 1 \
      -minMapQ 30 -minQ 20 \
      -out lava_gulls_thetas_group3_modern \
      -nThreads 12

thetaStat do_stat lava_gulls_thetas_group3_modern.thetas.idx

echo "Done! End time: $(date)"
