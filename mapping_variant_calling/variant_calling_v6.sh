#!/bin/bash

# Set working directory
WORK_DIR=~/Sequencing_Round1
cd $WORK_DIR

# Define input/output directories
BAM_DIR=${WORK_DIR}/complete_mapping
OUTPUT_DIR=${WORK_DIR}/angsd_variant_calling
mkdir -p $OUTPUT_DIR

# Define reference genome
REFERENCE=~/Doves_run/genome_data/yellow_legged_gull_genome.filtered.fna

# Define samples
SAMPLES=("LaGu_Workshop04" "LaGu_Workshop05" "LaGu_Workshop06" "LaGu_Workshop08" "LaGu_Workshop09")

# Define parameters for ANGSD
MIN_DEPTH=4         # Minimum depth parameter
MAX_DEPTH=40        # Maximum depth parameter
MIN_MAPQ=30         # Minimum mapping quality
MIN_QUAL=20         # Minimum base quality
MAF=0.1             # Minimum minor allele frequency
MAC=1               # Minimum minor allele count
THREADS=8           # Number of threads to use

# Create a list of BAM files for variant calling
BAMLIST=${OUTPUT_DIR}/bam_list.txt
> $BAMLIST  # Clear file if it exists

# Add all BAM files to the list - using dedup.bam files (with duplicates marked)
for SAMPLE in "${SAMPLES[@]}"; do
    if [ -f "${BAM_DIR}/${SAMPLE}.dedup.bam" ]; then
        echo "${BAM_DIR}/${SAMPLE}.dedup.bam" >> $BAMLIST
        echo "Added ${SAMPLE} to the analysis"
    else
        echo "WARNING: Dedup BAM file for ${SAMPLE} not found at ${BAM_DIR}/${SAMPLE}.dedup.bam"
        
        # Fallback to sorted BAM if dedup doesn't exist
        if [ -f "${BAM_DIR}/${SAMPLE}.sorted.bam" ]; then
            echo "${BAM_DIR}/${SAMPLE}.sorted.bam" >> $BAMLIST
            echo "Using sorted BAM for ${SAMPLE} instead"
        else
            echo "ERROR: No suitable BAM file found for ${SAMPLE}"
        fi
    fi
done

# Check if any BAM files were found
if [ ! -s "$BAMLIST" ]; then
    echo "ERROR: No valid BAM files found. Check the BAM_DIR path and sample names."
    exit 1
fi

echo "Starting ANGSD variant calling with parameters:"
echo "  Minimum mapping quality: $MIN_MAPQ"
echo "  Minimum base quality: $MIN_QUAL"
echo "  Minimum minor allele frequency: $MAF"
echo "  Number of samples: $(wc -l < $BAMLIST)"
echo "  BAM files list: $BAMLIST"

# Generate a BCF file based on suggested example from ANGSD output
echo "Generating variants with ANGSD using doBcf..."
angsd \
    -b $BAMLIST \
    -ref $REFERENCE \
    -out ${OUTPUT_DIR}/angsd_variants \
    -remove_bads 1 \
    -minMapQ $MIN_MAPQ \
    -minQ $MIN_QUAL \
    -minInd $(( $(wc -l < $BAMLIST) / 2 )) \
    -minMaf $MAF \
    -doBcf 1 \
    -doMajorMinor 1 \
    -doPost 1 \
    -gl 1 \
    -domaf 1 \
    -docounts 1 \
    -dogeno 1 \
    -P $THREADS

# Check if ANGSD ran successfully (now looking for BCF file)
if [ ! -f "${OUTPUT_DIR}/angsd_variants.bcf" ]; then
    echo "ERROR: ANGSD failed to generate BCF file."
    
    # Try one more time with absolute minimal parameters 
    echo "Retrying with absolute minimal parameters..."
    angsd \
        -b $BAMLIST \
        -ref $REFERENCE \
        -out ${OUTPUT_DIR}/angsd_minimal \
        -minMapQ $MIN_MAPQ \
        -minQ $MIN_QUAL \
        -GL 1 \
        -doMaf 1 \
        -doMajorMinor 1 \
        -SNP_pval 1e-6 \
        -P $THREADS
    
    echo "ANGSD completed. Generated basic statistics only."
    
    # Create a simple summary
    echo "Only basic statistics were generated." > ${OUTPUT_DIR}/summary.txt
    echo "Parameters used:" >> ${OUTPUT_DIR}/summary.txt
    echo "  Minimum mapping quality: $MIN_MAPQ" >> ${OUTPUT_DIR}/summary.txt
    echo "  Minimum base quality: $MIN_QUAL" >> ${OUTPUT_DIR}/summary.txt
    echo "  Number of samples: $(wc -l < $BAMLIST)" >> ${OUTPUT_DIR}/summary.txt
    echo "  BAM files used:" >> ${OUTPUT_DIR}/summary.txt
    cat $BAMLIST >> ${OUTPUT_DIR}/summary.txt
    
    echo "ANGSD completed with minimal variant statistics."
    echo "Output files are in: ${OUTPUT_DIR}"
    exit 0
fi

# Convert BCF to VCF and compress
echo "Converting BCF to VCF and compressing..."
bcftools view ${OUTPUT_DIR}/angsd_variants.bcf -Oz -o ${OUTPUT_DIR}/angsd_variants.vcf.gz

# Index the VCF file
echo "Indexing VCF file..."
bcftools index -t ${OUTPUT_DIR}/angsd_variants.vcf.gz

# Count number of variants in the VCF
echo "Counting variants..."
VAR_COUNT=$(bcftools view -H ${OUTPUT_DIR}/angsd_variants.vcf.gz | wc -l)
echo "Found $VAR_COUNT variants"

# Write summary file
echo "ANGSD variant calling complete!" > ${OUTPUT_DIR}/summary.txt
echo "Total variants called: $VAR_COUNT" >> ${OUTPUT_DIR}/summary.txt
echo "Parameters used:" >> ${OUTPUT_DIR}/summary.txt
echo "  Minimum mapping quality: $MIN_MAPQ" >> ${OUTPUT_DIR}/summary.txt
echo "  Minimum base quality: $MIN_QUAL" >> ${OUTPUT_DIR}/summary.txt
echo "  Minimum minor allele frequency: $MAF" >> ${OUTPUT_DIR}/summary.txt
echo "  Number of samples: $(wc -l < $BAMLIST)" >> ${OUTPUT_DIR}/summary.txt
echo "  BAM files used:" >> ${OUTPUT_DIR}/summary.txt
cat $BAMLIST >> ${OUTPUT_DIR}/summary.txt

echo "ANGSD variant calling completed successfully!"
echo "VCF file is available at: ${OUTPUT_DIR}/angsd_variants.vcf.gz"
