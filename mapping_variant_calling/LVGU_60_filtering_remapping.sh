#!/bin/bash

# Process new LeucoFulig ONT data and map to ONT-filtered Yellow-legged Gull reference

LONGREAD_FASTQ=~/Sequencing_Mapping/ref_genome_LVGU60/LeucoFulig_basecall_pass.fastq.gz
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna
OUTPUT_DIR=~/Sequencing_Combined/modern_remap_paleomix_Oct25

mkdir -p $OUTPUT_DIR

# 1. Adapter trimming
porechop -i $LONGREAD_FASTQ -o temp_trimmed.fastq.gz --threads 8 --discard_middle

# 2. Quality filtering
gunzip -c temp_trimmed.fastq.gz | NanoFilt -q 18 -l 2000 --headcrop 50 --tailcrop 50 | gzip > temp_filtered.fastq.gz

# 3. Map and filter to ONT-filtered reference
minimap2 -ax map-ont --secondary=no -N 1 --min-chain-score 1000 \
    -R '@RG\tID:LVGU_60_fastq\tSM:LVGU_60_fastq\tLB:lib1\tPL:ONT' \
    $REFERENCE temp_filtered.fastq.gz \
    | samtools view -q 20 -F 4 -b \
    | samtools sort -o $OUTPUT_DIR/LVGU_60_fastq.YellowLeggedGull_ONT.bam

# 4. Index and cleanup
samtools index $OUTPUT_DIR/LVGU_60_fastq.YellowLeggedGull_ONT.bam
rm temp_*.fastq.gz

echo "Done. Ready for ANGSD with ONT-filtered reference."
samtools flagstat $OUTPUT_DIR/LVGU_60_fastq.YellowLeggedGull_ONT.bam
