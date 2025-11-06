#!/bin/bash

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
BAM_LIST=~/modern_samples_only_bam_list.txt
REF=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna

cd $WORK_DIR

echo "=========================================="
echo "Group 4a: Scaffolds (modern samples, n=33)"
echo "=========================================="
echo "Start time: $(date)"

echo "Step 1: Generating SAF..."
angsd -bam $BAM_LIST \
      -ref $REF \
      -anc $REF \
      -doSaf 1 -GL 1 \
      -minMapQ 30 -minQ 20 \
      -rf ~/regions_group4a.txt \
      -out lava_gulls_saf_group4a_modern \
      -nThreads 12

echo "SAF complete!"

echo "Step 2: Calculating SFS..."
realSFS lava_gulls_saf_group4a_modern.saf.idx > lava_gulls_group4a_modern.sfs

echo "SFS complete!"

echo "Step 3: Calculating thetas..."
realSFS saf2theta lava_gulls_saf_group4a_modern.saf.idx \
                   -sfs lava_gulls_group4a_modern.sfs \
                   -outname lava_gulls_thetas_group4a_modern

echo "Thetas complete!"

echo "Step 4: Summary statistics..."
thetaStat do_stat lava_gulls_thetas_group4a_modern.thetas.idx

echo "=========================================="
echo "GROUP 4a COMPLETE! End time: $(date)"
echo "=========================================="
