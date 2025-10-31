#!/bin/bash

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/all_samples_bam_list_v2.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

cd $WORK_DIR

echo "Calculating thetas (π, Watterson's theta, Tajima's D)..."
echo "Start time: $(date)"

angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -doThetas 1 -doSaf 1 \
      -pest lava_gulls_all_merged.sfs \
      -GL 1 \
      -minMapQ 30 -minQ 20 \
      -out lava_gulls_all_merged_thetas \
      -nThreads 12

echo "Thetas calculated!"

echo "Calculating summary statistics..."
thetaStat do_stat lava_gulls_all_merged_thetas.thetas.idx

echo "Calculating windowed statistics (10kb windows, 5kb step)..."
thetaStat do_stat lava_gulls_all_merged_thetas.thetas.idx \
          -win 10000 -step 5000 \
          -outnames lava_gulls_all_merged_thetas_10kb.thetasWindow

echo "ALL DONE!"
echo "End time: $(date)"
echo ""
echo "Results files:"
echo "  - lava_gulls_all_merged_thetas.thetas.idx.pestPG (genome-wide statistics)"
echo "  - lava_gulls_all_merged_thetas_10kb.thetasWindow.pestPG (windowed statistics)"
