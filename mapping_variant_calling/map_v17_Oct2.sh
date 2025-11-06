#!/bin/bash

REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna
R1=~/Sequencing_Combined/cleaned_fastq/LVGU_10_R1_cleaned.fq.gz
R2=~/Sequencing_Combined/cleaned_fastq/LVGU_10_R2_cleaned.fq.gz

minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
samtools view -bS -q 20 -@ 8 - | \
samtools sort -@ 8 -o ~/Sequencing_Combined/complete_mapping/LVGU_10_testOct_temp.bam

samtools markdup -r -@ 8 \
  ~/Sequencing_Combined/complete_mapping/LVGU_10_testOct_temp.bam \
  ~/Sequencing_Combined/complete_mapping/LVGU_10_testOct.bam

samtools index ~/Sequencing_Combined/complete_mapping/LVGU_10_testOct.bam
rm ~/Sequencing_Combined/complete_mapping/LVGU_10_testOct_temp.bam
