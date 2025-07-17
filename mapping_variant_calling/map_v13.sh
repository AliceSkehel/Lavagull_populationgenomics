#!/bin/bash

# Debug mapping script with explicit checkpoints for MC tag troubleshooting
READS_DIR=~/Sequencing_Round1/fastp_minimal_cleaned
REF_GENOME=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna
OUTPUT_DIR=~/Sequencing_Round1/complete_mapping
DEBUG_DIR=$OUTPUT_DIR/debug

mkdir -p $DEBUG_DIR

# Process one sample for debugging
SAMPLE="LaGu_Workshop04"  # Change this to test different samples
R1=$READS_DIR/${SAMPLE}_R1_minimal.fq.gz
R2=$READS_DIR/${SAMPLE}_R2_minimal.fq.gz

echo "=== DEBUG PROCESSING: $SAMPLE ==="
echo "Input files:"
echo "  R1: $R1"
echo "  R2: $R2"

# Define intermediate files
INITIAL_SAM=$DEBUG_DIR/${SAMPLE}.01_initial.sam
FILTERED_BAM=$DEBUG_DIR/${SAMPLE}.02_filtered.bam
NAMESORTED_BAM=$DEBUG_DIR/${SAMPLE}.03_namesorted.bam
FIXMATE_BAM=$DEBUG_DIR/${SAMPLE}.04_fixmate.bam
COORDSORTED_BAM=$DEBUG_DIR/${SAMPLE}.05_coordsorted.bam
FINAL_BAM=$DEBUG_DIR/${SAMPLE}.06_final.bam

echo ""
echo "=== STEP 1: ALIGNMENT ==="
echo "Command: minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 > $INITIAL_SAM"
minimap2 -ax sr -t 8 $REF_GENOME $R1 $R2 > $INITIAL_SAM

if [ $? -eq 0 ]; then
    echo "✅ STEP 1 SUCCESS: Alignment completed"
    echo "   Output: $INITIAL_SAM"
    # Quick stats
    echo "   Total alignments: $(wc -l < $INITIAL_SAM)"
else
    echo "❌ STEP 1 FAILED: Alignment failed"
    exit 1
fi

echo ""
echo "=== STEP 2: FILTERING (MAPQ >= 20) ==="
echo "Command: samtools view -bS -q 20 $INITIAL_SAM -o $FILTERED_BAM"
samtools view -bS -q 20 $INITIAL_SAM -o $FILTERED_BAM

if [ $? -eq 0 ]; then
    echo "✅ STEP 2 SUCCESS: Quality filtering completed"
    echo "   Output: $FILTERED_BAM"
    # Check file integrity
    samtools quickcheck $FILTERED_BAM && echo "   File integrity: OK" || echo "   File integrity: FAILED"
    # Basic stats
    echo "   Reads after filtering: $(samtools view -c $FILTERED_BAM)"
else
    echo "❌ STEP 2 FAILED: Quality filtering failed"
    exit 1
fi

echo ""
echo "=== STEP 3: NAME SORTING ==="
echo "Command: samtools sort -n $FILTERED_BAM -o $NAMESORTED_BAM"
samtools sort -n $FILTERED_BAM -o $NAMESORTED_BAM

if [ $? -eq 0 ]; then
    echo "✅ STEP 3 SUCCESS: Name sorting completed"
    echo "   Output: $NAMESORTED_BAM"
    # Verify sort order
    samtools view -H $NAMESORTED_BAM | grep "SO:" || echo "   No sort order in header"
    echo "   First few read names:"
    samtools view $NAMESORTED_BAM | head -3 | cut -f1
else
    echo "❌ STEP 3 FAILED: Name sorting failed"
    exit 1
fi

echo ""
echo "=== STEP 4: FIXMATE (CRITICAL FOR MC TAGS) ==="
echo "Command: samtools fixmate -m $NAMESORTED_BAM $FIXMATE_BAM"
samtools fixmate -m $NAMESORTED_BAM $FIXMATE_BAM

if [ $? -eq 0 ]; then
    echo "✅ STEP 4 SUCCESS: Fixmate completed"
    echo "   Output: $FIXMATE_BAM"
    # CRITICAL CHECK: Verify MC tags were added
    MC_COUNT=$(samtools view $FIXMATE_BAM | head -1000 | grep -c "MC:Z:")
    MS_COUNT=$(samtools view $FIXMATE_BAM | head -1000 | grep -c "ms:i:")
    echo "   MC tags found in first 1000 reads: $MC_COUNT"
    echo "   ms tags found in first 1000 reads: $MS_COUNT"
    
    if [ $MC_COUNT -gt 0 ]; then
        echo "   ✅ MC tags successfully added!"
        # Show example MC tag
        echo "   Example MC tag: $(samtools view $FIXMATE_BAM | grep "MC:Z:" | head -1 | grep -o "MC:Z:[^[:space:]]*")"
    else
        echo "   ❌ WARNING: No MC tags found! This will cause markdup to fail."
    fi
else
    echo "❌ STEP 4 FAILED: Fixmate failed"
    exit 1
fi

echo ""
echo "=== STEP 5: COORDINATE SORTING ==="
echo "Command: samtools sort $FIXMATE_BAM -o $COORDSORTED_BAM"
samtools sort $FIXMATE_BAM -o $COORDSORTED_BAM

if [ $? -eq 0 ]; then
    echo "✅ STEP 5 SUCCESS: Coordinate sorting completed"
    echo "   Output: $COORDSORTED_BAM"
    # Verify sort order
    samtools view -H $COORDSORTED_BAM | grep "SO:" || echo "   No sort order in header"
    # Verify MC tags survived sorting
    MC_COUNT_AFTER=$(samtools view $COORDSORTED_BAM | head -1000 | grep -c "MC:Z:")
    echo "   MC tags after sorting: $MC_COUNT_AFTER"
else
    echo "❌ STEP 5 FAILED: Coordinate sorting failed"
    exit 1
fi

echo ""
echo "=== STEP 6: DUPLICATE MARKING ==="
echo "Command: samtools markdup -r $COORDSORTED_BAM $FINAL_BAM"
samtools markdup -r $COORDSORTED_BAM $FINAL_BAM

if [ $? -eq 0 ]; then
    echo "✅ STEP 6 SUCCESS: Duplicate marking completed"
    echo "   Output: $FINAL_BAM"
    # Final stats
    echo "   Final read count: $(samtools view -c $FINAL_BAM)"
    echo "   Creating index..."
    samtools index $FINAL_BAM
    echo "   ✅ Index created: ${FINAL_BAM}.bai"
else
    echo "❌ STEP 6 FAILED: Duplicate marking failed"
    echo "   This is likely due to missing MC tags from Step 4"
    exit 1
fi

echo ""
echo "=== FINAL SUMMARY ==="
echo "✅ All steps completed successfully!"
echo "Debug files saved in: $DEBUG_DIR"
echo "Final output: $FINAL_BAM"
echo ""
echo "Intermediate files created:"
echo "  1. $INITIAL_SAM (initial alignment)"
echo "  2. $FILTERED_BAM (quality filtered)"
echo "  3. $NAMESORTED_BAM (name sorted for fixmate)"
echo "  4. $FIXMATE_BAM (MC tags added - CHECK THIS ONE!)"
echo "  5. $COORDSORTED_BAM (coordinate sorted for markdup)"
echo "  6. $FINAL_BAM (final deduplicated file)"
echo ""
echo "To clean up debug files: rm $DEBUG_DIR/${SAMPLE}.*"
