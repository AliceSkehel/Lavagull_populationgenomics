#!/bin/bash

# Simple fast ONT mapping
FASTQ_DIR="/home/askehel/mitogenome/mito_fastq_merged"
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna"
OUTPUT_DIR="/home/askehel/mitogenome/mapping_results"

mkdir -p ${OUTPUT_DIR}

for fastq in ${FASTQ_DIR}/*.fastq.gz; do
    name=$(basename ${fastq} .fastq.gz)
    echo "Processing: ${name}"
    
    (
    minimap2 -ax map-ont -t 3 ${REFERENCE} ${fastq} | \
    samtools view -bS -q 20 -F 256 -@ 3 - | \
    samtools sort -@ 3 - | \
    samtools markdup -r -@ 3 - ${OUTPUT_DIR}/${name}.bam
    
    samtools index -@ 3 ${OUTPUT_DIR}/${name}.bam
    ) &
    
    # Limit to 4 parallel jobs
    if (( $(jobs -r | wc -l) >= 4 )); then
        wait -n
    fi
done

wait
echo "Done!"

