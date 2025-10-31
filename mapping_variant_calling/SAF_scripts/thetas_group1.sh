#!/bin/bash

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/all_samples_bam_list_v2.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

cd $WORK_DIR

echo "Calculating thetas for Group 1..."
echo "Start time: $(date)"

# Calculate SFS for group 1
realSFS lava_gulls_saf_group1.saf.idx > lava_gulls_group1.sfs

# Calculate thetas
angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -rf ~/regions_group1.txt \
      -doThetas 1 -doSaf 1 \
      -pest lava_gulls_group1.sfs \
      -GL 1 \
      -minMapQ 30 -minQ 20 \
      -out lava_gulls_thetas_group1 \
      -nThreads 12

echo "Calculating summary statistics..."
thetaStat do_stat lava_gulls_thetas_group1.thetas.idx

echo "Done! End time: $(date)"
