#!/bin/bash

# BWA mapping script for Lava Gull mtDNA
# Usage: ./map_lava_gull_mtdna.sh

set -e  # Exit on any error

REFERENCE_DIR="/home/askehel/mitogenome"
REFERENCE_FILE="$REFERENCE_DIR/KM507782_mtDNA.fasta"
DATA_DIR="$REFERENCE_DIR/LavaGulls_mtDNA"
OUTPUT_DIR="$REFERENCE_DIR/bwa_mappings"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if reference exists
if [ ! -f "$REFERENCE_FILE" ]; then
    echo "Error: Reference file not found at $REFERENCE_FILE"
    echo "Please download the reference first"
    exit 1
fi

echo "=== BWA Mapping Pipeline for Lava Gull mtDNA ==="
echo "Reference: $REFERENCE_FILE"
echo "Data directory: $DATA_DIR"
echo "Output directory: $OUTPUT_DIR"
echo

# Step 1: Index the reference genome
echo "Step 1: Indexing reference genome..."
bwa index "$REFERENCE_FILE"
echo "Reference indexing complete."
echo

# Step 2: Map each FASTA file
echo "Step 2: Mapping sequences..."

for fasta_file in "$DATA_DIR"/*.fa; do
    if [ -f "$fasta_file" ]; then
        # Extract base name without path and extension
        base_name=$(basename "$fasta_file" .fa)
        
        echo "Processing: $base_name"
        
        # BWA-MEM alignment (good for longer sequences like complete mtDNA)
        bwa mem -t 4 "$REFERENCE_FILE" "$fasta_file" > "$OUTPUT_DIR/${base_name}_vs_ref.sam"
        
        # Convert SAM to BAM
        samtools view -bS "$OUTPUT_DIR/${base_name}_vs_ref.sam" > "$OUTPUT_DIR/${base_name}_vs_ref.bam"
        
        # Sort BAM file
        samtools sort "$OUTPUT_DIR/${base_name}_vs_ref.bam" -o "$OUTPUT_DIR/${base_name}_vs_ref_sorted.bam"
        
        # Index sorted BAM
        samtools index "$OUTPUT_DIR/${base_name}_vs_ref_sorted.bam"
        
        # Generate mapping statistics
        samtools flagstat "$OUTPUT_DIR/${base_name}_vs_ref_sorted.bam" > "$OUTPUT_DIR/${base_name}_mapping_stats.txt"
        
        # Clean up intermediate files
        rm "$OUTPUT_DIR/${base_name}_vs_ref.sam" "$OUTPUT_DIR/${base_name}_vs_ref.bam"
        
        echo "  -> Completed: ${base_name}_vs_ref_sorted.bam"
    fi
done

echo
echo "=== Mapping Summary ==="
echo "Results saved in: $OUTPUT_DIR"
echo "Files generated per sample:"
echo "  - *_vs_ref_sorted.bam (sorted alignment)"
echo "  - *_vs_ref_sorted.bam.bai (index)"
echo "  - *_mapping_stats.txt (alignment statistics)"
echo
echo "To view mapping stats:"
echo "cat $OUTPUT_DIR/*_mapping_stats.txt"
