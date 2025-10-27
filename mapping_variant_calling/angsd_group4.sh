#!/bin/bash
# ANGSD Group 4: Chr11+ (242 remaining chromosomes, ~200 Mb)

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35
mkdir -p $WORK_DIR

echo "Starting ANGSD Group 4 (Chr11+)..."
echo "Start time: $(date)"

angsd -bam ~/all_samples_bam_list_v2.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_all35_group4 \
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
      -rf ~/regions_group4.txt \
      -SNP_pval 1e-6 \
      -nThreads 12

echo "Done! End time: $(date)"
echo "SNP count:"
zcat $WORK_DIR/lava_gulls_all35_group4.mafs.gz | wc -l
