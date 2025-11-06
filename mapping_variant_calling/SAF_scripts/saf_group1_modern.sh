#!/bin/bash

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/modern_samples_only_bam_list.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

mkdir -p $WORK_DIR
cd $WORK_DIR

echo "=========================================="
echo "Group 1: Chr1 (modern samples, n=33)"
echo "=========================================="
echo "Start time: $(date)"

# Step 1: Generate SAF
echo "Step 1: Generating SAF..."
angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -doSaf 1 -GL 1 \
      -minMapQ 30 -minQ 20 \
      -rf ~/regions_group1.txt \
      -out lava_gulls_saf_group1_modern \
      -nThreads 12

echo "SAF complete!"

# Step 2: Calculate SFS
echo "Step 2: Calculating SFS..."
realSFS lava_gulls_saf_group1_modern.saf.idx > lava_gulls_group1_modern.sfs

echo "SFS complete!"

# Step 3: Calculate Thetas
echo "Step 3: Calculating thetas..."
angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -rf ~/regions_group1.txt \
      -doThetas 1 -doSaf 1 \
      -pest lava_gulls_group1_modern.sfs \
      -GL 1 \
      -minMapQ 30 -minQ 20 \
      -out lava_gulls_thetas_group1_modern \
      -nThreads 12

echo "Thetas complete!"

# Step 4: Summary statistics
echo "Step 4: Summary statistics..."
thetaStat do_stat lava_gulls_thetas_group1_modern.thetas.idx

echo "=========================================="
echo "GROUP 1 COMPLETE!"
echo "End time: $(date)"
echo "=========================================="
