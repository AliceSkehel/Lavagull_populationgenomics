#!/bin/bash

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/modern_samples_only_bam_list.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

cd $WORK_DIR

echo "Calculating thetas for Group 4c (modern samples only, n=33)..."
echo "Start time: $(date)"

realSFS lava_gulls_saf_group4c.saf.idx > lava_gulls_group4c_modern.sfs

angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -rf ~/regions_group4c.txt \
      -doThetas 1 -doSaf 1 \
      -pest lava_gulls_group4c_modern.sfs \
      -GL 1 \
      -minMapQ 30 -minQ 20 \
      -out lava_gulls_thetas_group4c_modern \
      -nThreads 12

thetaStat do_stat lava_gulls_thetas_group4c_modern.thetas.idx

echo "Done! End time: $(date)"
