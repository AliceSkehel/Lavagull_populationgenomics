#!/bin/bash
# SAF Group 1: Chr1 only (224 Mb)

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/all_samples_bam_list_v2.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

mkdir -p $WORK_DIR
cd $WORK_DIR

echo "Starting SAF Group 1 (Chr1)..."
echo "Start time: $(date)"

angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -doSaf 1 -GL 1 \
      -minMapQ 30 -minQ 20 \
      -rf ~/regions_group1.txt \
      -out lava_gulls_saf_group1 \
      -nThreads 12

echo "Done! End time: $(date)"
