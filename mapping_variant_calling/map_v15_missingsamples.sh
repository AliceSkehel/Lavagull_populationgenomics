#!/bin/bash

# Map only missing samples - no parallel processing
READS_DIR=~/Sequencing_Round2/fastp_polyg_cleaned/cleaned_fastq
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.fna
OUTPUT_DIR=~/Sequencing_Round2/complete_mapping/unfiltered_ref

mkdir -p $OUTPUT_DIR

echo "Mapping only missing samples (sequential processing)..."
echo "======================================================="

for R1_FILE in $READS_DIR/*_R1_cleaned.fq.gz; do
    SAMPLE=$(basename $R1_FILE _R1_cleaned.fq.gz)
    
    # Skip if already done
    if [ -f "$OUTPUT_DIR/${SAMPLE}.bam" ] && [ -f "$OUTPUT_DIR/${SAMPLE}.bam.bai" ]; then
        echo "✅ $SAMPLE already completed, skipping..."
        continue
    fi
    
    echo ""
    echo "🔄 Processing: $SAMPLE"
    echo "Time: $(date)"
    
    R1=$READS_DIR/${SAMPLE}_R1_cleaned.fq.gz
    R2=$READS_DIR/${SAMPLE}_R2_cleaned.fq.gz
    
    # Check input files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "❌ Missing input files for $SAMPLE"
        continue
    fi
    
    echo "  Input: $(basename $R1) + $(basename $R2)"
    
    # Map and process with error checking
    echo "  🗺️  Mapping..."
    if minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 | \
       samtools view -bS -q 20 -@ 8 - | \
       samtools sort -n -@ 8 - | \
       samtools fixmate -m -@ 8 - - | \
       samtools sort -@ 8 - | \
       samtools view -h -@ 8 - | \
       awk '/^@/ || /MC:Z:/' | \
       samtools view -b -@ 8 - | \
       samtools markdup -r -@ 8 - $OUTPUT_DIR/${SAMPLE}.bam; then
        
        echo "  📊 Indexing..."
        samtools index -@ 8 $OUTPUT_DIR/${SAMPLE}.bam
        
        echo "  ✅ Success: ${SAMPLE}.bam"
        
        # Quick stats
        echo "  📈 Quick stats:"
        samtools flagstat $OUTPUT_DIR/${SAMPLE}.bam | head -3
        
    else
        echo "  ❌ Failed: $SAMPLE"
        # Remove failed BAM if it exists
        rm -f $OUTPUT_DIR/${SAMPLE}.bam
    fi
done

echo ""
echo "======================================================="
echo "🎉 Mapping complete!"
echo ""
echo "Final count:"
BAM_COUNT=$(ls $OUTPUT_DIR/*.bam 2>/dev/null | wc -l)
echo "Total BAM files: $BAM_COUNT"
echo ""
echo "All samples:"
ls -lh $OUTPUT_DIR/*.bam 2>/dev/null | tail -10
