#!/bin/bash

# ANGSD v8 - Illumina Only with SNP Filtering (32 samples)
# Adds -SNP_pval filter to get high-confidence SNPs (matches killer whale paper approach)

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling
mkdir -p $WORK_DIR

echo "Starting ANGSD variant calling with SNP filters..."
echo "Start time: $(date)"

# ANGSD with illumina-only samples + SNP filtering
angsd -bam ~/illumina_only_bam_list.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_illumina_v8 \
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
      -minInd 25 \
      -doCounts 1 \
      -setMinDepth 120 \
      -setMaxDepth 2534 \
      -setMinDepthInd 3 \
      -setMaxDepthInd 30 \
      -rf ~/regions_no_mtdna.txt \
      -SNP_pval 1e-6 \
      -nThreads 12

echo "Done!"
echo "End time: $(date)"

# Print summary
echo ""
echo "Output files created:"
ls -lh $WORK_DIR/lava_gulls_illumina_v8.*

echo ""
echo "To check number of SNPs:"
echo "zcat $WORK_DIR/lava_gulls_illumina_v8.mafs.gz | wc -l"

