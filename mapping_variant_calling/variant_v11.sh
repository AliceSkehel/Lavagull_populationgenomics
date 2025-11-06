#!/bin/bash

# ANGSD v11 - All 35 samples: 32 modern + 2 ancient + 1 LVGU_60_fastq
# Adjusted for low-coverage ancient samples (1.3X, 2.3X) using genotype likelihoods
# Depth calculations: min = 35 × 1.3 × 0.5 = 23, max = 35 × 29.1 × 3 = 3060

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35
mkdir -p $WORK_DIR

echo "Starting ANGSD variant calling with all 35 samples (including ancient)..."
echo "Start time: $(date)"

# ANGSD with genotype likelihoods for mixed coverage (1.3X to 29.1X)
angsd -bam ~/all_samples_bam_list_v2.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_all35_v11 \
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
      -C 50 \
      -minInd 28 \
      -doCounts 1 \
      -setMinDepth 23 \
      -setMaxDepth 3060 \
      -setMinDepthInd 1 \
      -setMaxDepthInd 60 \
      -rf ~/regions_no_mtdna.txt \
      -SNP_pval 1e-6 \
      -nThreads 12

echo "Done!"
echo "End time: $(date)"

echo ""
echo "Output files created:"
ls -lh $WORK_DIR/lava_gulls_all35_v11.*

echo ""
echo "To check number of SNPs:"
echo "zcat $WORK_DIR/lava_gulls_all35_v11.mafs.gz | wc -l"
