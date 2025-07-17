!/bin/bash
# Fixed Illumina mapping script with proper paired-end processing - map_v10.sh
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping
TEMP_DIR=~/tmp
mkdir -p $OUTPUT_DIR
mkdir -p $TEMP_DIR
export TMPDIR=$TEMP_DIR

for R1_FILE in $READS_DIR/*_R1_minimal.fq.gz; do
    SAMPLE=$(basename $R1_FILE _R1_minimal.fq.gz)
    R1=$READS_DIR/${SAMPLE}_R1_minimal.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_minimal.fq.gz
    
    echo "Processing: $SAMPLE"
    
    (
    # Fixed pipeline with proper paired-end processing
    minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
    samtools view -bS -q 20 -@ 8 - | \
    samtools sort -n -@ 8 -T $TEMP_DIR/${SAMPLE}_sort1 - | \
    samtools fixmate -m -@ 8 - - | \
    samtools sort -@ 8 -T $TEMP_DIR/${SAMPLE}_sort2 - | \
    samtools markdup -r -@ 8 - $OUTPUT_DIR/${SAMPLE}.bam
    
    samtools index -@ 8 $OUTPUT_DIR/${SAMPLE}.bam
    echo "Completed: $SAMPLE"
    ) &
done

wait
echo "All samples processed!"
echo "Output files are in: $OUTPUT_DIR"
