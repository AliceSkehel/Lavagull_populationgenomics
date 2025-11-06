#!/bin/bash

# ANGSD v10 - 33 modern samples only (excluding ancient LVGU_473, LVGU_476)

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_modern33
mkdir -p $WORK_DIR

echo "Starting ANGSD variant calling with 33 modern samples..."
echo "Start time: $(date)"

angsd -bam ~/paleomix_bam_list_33modern.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_modern33 \
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
      -setMinDepth 129 \
      -setMaxDepth 2882 \
      -setMinDepthInd 3 \
      -setMaxDepthInd 87 \
      -rf ~/regions_no_mtdna.txt \
      -SNP_pval 1e-6 \
      -nThreads 12

echo "Done!"
echo "End time: $(date)"

echo ""
echo "Output files created:"
ls -lh $WORK_DIR/lava_gulls_modern33.*

echo ""
echo "To check number of SNPs:"
echo "zcat $WORK_DIR/lava_gulls_modern33.mafs.gz | wc -l"
