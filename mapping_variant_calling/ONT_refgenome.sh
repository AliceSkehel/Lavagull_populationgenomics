#!/bin/bash


# Stringent long-read processing for pop gen
LONGREAD_FASTQ=~/Sequencing_Round2/fastq_ref_genome/LavaGull_reference_sample.fastq.gz
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.fna
OUTPUT_DIR=~/Sequencing_Round2/complete_mapping

# 1. Adapter trimming
porechop -i $LONGREAD_FASTQ -o temp_trimmed.fastq.gz --threads 8 --discard_middle

# 2. Quality filtering (decompress first)
gunzip -c temp_trimmed.fastq.gz | NanoFilt -q 18 -l 2000 --headcrop 50 --tailcrop 50 > temp_filtered.fastq.gz

# 3. Map and filter
minimap2 -ax map-ont --secondary=no -N 1 --min-chain-score 1000 \
    -R '@RG\tID:LVGU_longread\tSM:LVGU_longread\tLB:lib1\tPL:ONT' \
    $REFERENCE temp_filtered.fastq.gz \
    | samtools view -q 20 -F 4 -b \
    | samtools sort -o $OUTPUT_DIR/LVGU_longread.bam

# 4. Index and cleanup
samtools index $OUTPUT_DIR/LVGU_longread.bam
rm temp_*.fastq.gz

echo "Done. Ready for ANGSD."
samtools flagstat $OUTPUT_DIR/LVGU_longread.bam
