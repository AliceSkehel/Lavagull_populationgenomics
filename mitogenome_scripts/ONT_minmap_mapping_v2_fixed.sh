#!/bin/bash
# Reprocess failed ONT samples with temp directory fix

FASTQ_DIR="/home/askehel/mitogenome/mito_fastq_merged"
REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna"
OUTPUT_DIR="/home/askehel/mitogenome/mapping_results"
TEMP_DIR="/home/askehel/tmp"

# Create directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}
export TMPDIR=${TEMP_DIR}

# Define the failed samples
FAILED_SAMPLES=(
    "CF_LG_barcode11_merged_porechop-abi"
    "LG_GP_barcode01_merged_porechop-abi"
    "LG_GP_barcode02_merged_porechop-abi" 
    "LG_GP_barcode03_merged_porechop-abi"
)

echo "Removing old failed BAM files..."
for sample in "${FAILED_SAMPLES[@]}"; do
    if [ -f "${OUTPUT_DIR}/${sample}.bam" ]; then
        rm "${OUTPUT_DIR}/${sample}.bam"
        echo "Removed: ${sample}.bam"
    fi
    if [ -f "${OUTPUT_DIR}/${sample}.bam.bai" ]; then
        rm "${OUTPUT_DIR}/${sample}.bam.bai"
        echo "Removed: ${sample}.bam.bai"
    fi
done

echo "Starting reprocessing of failed samples..."

for sample in "${FAILED_SAMPLES[@]}"; do
    fastq="${FASTQ_DIR}/${sample}.fastq.gz"
    
    if [ -f "$fastq" ]; then
        echo "Processing: ${sample}"
        
        (
        minimap2 -ax map-ont -t 3 ${REFERENCE} ${fastq} | \
        samtools view -bS -q 20 -F 256 -@ 3 - | \
        samtools sort -@ 3 -T ${TEMP_DIR}/${sample}_sort - | \
        samtools markdup -r -@ 3 - ${OUTPUT_DIR}/${sample}.bam
        
        samtools index -@ 3 ${OUTPUT_DIR}/${sample}.bam
        echo "Completed: ${sample}"
        ) &
        
        # Limit to 2 parallel jobs to reduce resource competition
        if (( $(jobs -r | wc -l) >= 2 )); then
            wait -n
        fi
    else
        echo "Warning: FASTQ file not found: $fastq"
    fi
done

wait
echo "Reprocessing complete!"
echo "Checking results..."

# Verify the results
for sample in "${FAILED_SAMPLES[@]}"; do
    bam_file="${OUTPUT_DIR}/${sample}.bam"
    if [ -f "$bam_file" ]; then
        size=$(ls -lh "$bam_file" | awk '{print $5}')
        echo "${sample}.bam: $size"
    else
        echo "${sample}.bam: NOT FOUND"
    fi
done
