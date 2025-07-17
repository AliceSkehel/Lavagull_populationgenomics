#!/bin/bash

# Map FASTQ files to reference genome
FASTQ_DIR="/home/askehel/mitogenome/mito_fastq_merged"
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna"
OUTPUT_DIR="/home/askehel/mitogenome/mapping_results"

mkdir -p ${OUTPUT_DIR}

for fastq in ${FASTQ_DIR}/*.fastq.gz; do
    name=$(basename ${fastq} .fastq.gz)
    echo "Mapping: ${name}"
    
    minimap2 -ax map-ont ${REFERENCE} ${fastq} | \
    samtools view -bS -q 20 -F 256 - | \
    samtools sort - | \
    samtools markdup -r - ${OUTPUT_DIR}/${name}.bam
    
    samtools index ${OUTPUT_DIR}/${name}.bam
done

echo "Done! BAM files in: ${OUTPUT_DIR}"



