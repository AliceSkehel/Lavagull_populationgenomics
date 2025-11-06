#!/bin/bash

# ANGSD v9 - All samples: 32 modern + 2 ancient + 1 LVGU_60 (35 total)
# SNP filtering with -SNP_pval 1e-6 only

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling
mkdir -p $WORK_DIR

echo "Starting ANGSD variant calling with all 35 samples..."
echo "Start time: $(date)"

# ANGSD with all samples + SNP filtering
angsd -bam ~/all_samples_bam_list.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_all35 \
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
      -setMinDepth 140 \
      -setMaxDepth 2800 \
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
ls -lh $WORK_DIR/lava_gulls_all35.*

echo ""
echo "To check number of SNPs:"
echo "zcat $WORK_DIR/lava_gulls_all35.mafs.gz | wc -l"
