
#!/bin/bash

# Re-map ONT data to match Illumina reference genome
# This fixes the reference genome mismatch causing alignment issues

REFERENCE="/home/askehel/Doves_run/genome_data/yellow_legged_gull_genome.ont_filtered.fna"
FASTQ_DIR="/home/askehel/mitogenome/mito_fastq_merged"
OUTPUT_DIR="/home/askehel/mitogenome/corrected_mapping"

mkdir -p ${OUTPUT_DIR}

echo "Re-mapping ONT samples to correct reference genome..."
echo "Reference: ${REFERENCE}"

# Re-map all ONT samples to the correct reference
for fastq in ${FASTQ_DIR}/LVGU_*.fastq.gz; do
    if [ -f "$fastq" ]; then
        sample=$(basename $fastq .fastq.gz)
        echo "Processing: $sample"
        
        minimap2 -ax map-ont -t 8 ${REFERENCE} ${fastq} | \
        samtools sort -@ 8 -o ${OUTPUT_DIR}/${sample}_ONT_corrected.bam
        
        samtools index ${OUTPUT_DIR}/${sample}_ONT_corrected.bam
        
        # Check output size
        size=$(ls -lh ${OUTPUT_DIR}/${sample}_ONT_corrected.bam | awk '{print $5}')
        echo "  Completed: ${sample}_ONT_corrected.bam (${size})"
    fi
done

echo "Re-mapping complete!"
echo "Output directory: ${OUTPUT_DIR}"
