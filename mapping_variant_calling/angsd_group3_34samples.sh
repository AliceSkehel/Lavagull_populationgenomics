#!/bin/bash
# ANGSD Group 3 - 34 samples

WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_34samples
mkdir -p $WORK_DIR

echo "Starting ANGSD Group 3 - 34 samples..."
echo "Start time: $(date)"

angsd -bam ~/all_samples_bam_list_no_LVGU60.txt \
      -ref ~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna \
      -out $WORK_DIR/lava_gulls_34samples_group3 \
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
      -minInd 27 \
      -doCounts 1 \
      -setMinDepth 22 \
      -setMaxDepth 2960 \
      -setMinDepthInd 1 \
      -setMaxDepthInd 60 \
      -rf ~/regions_group3.txt \
      -SNP_pval 1e-6 \
      -nThreads 12

echo "Done! End time: $(date)"
echo "SNP count:"
zcat $WORK_DIR/lava_gulls_34samples_group3.mafs.gz | wc -l
