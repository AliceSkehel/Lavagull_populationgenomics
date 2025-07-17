#!/bin/bash

# Navigate to your directory
cd mitogenome/mito_fastq_merged/

# Create backup directory (optional but recommended)
mkdir -p ../backup_merged
cp *.fastq.gz ../backup_merged/

# Rename files using mv command
mv "CF_LG_barcode07_merged_porechop-abi.fastq.gz" "LVGU_030.fastq.gz"
mv "CF_LG_barcode08_merged_porechop-abi.fastq.gz" "LVGU_034.fastq.gz"
mv "CF_LG_barcode09_merged_porechop-abi.fastq.gz" "LVGU_038.fastq.gz"
mv "CF_LG_barcode10_merged_porechop-abi.fastq.gz" "LVGU_039.fastq.gz"
mv "CF_LG_barcode11_merged_porechop-abi.fastq.gz" "LVGU_042.fastq.gz"
mv "LG_GP_barcode01_merged_porechop-abi.fastq.gz" "LVGU_041.fastq.gz"
mv "LG_GP_barcode02_merged_porechop-abi.fastq.gz" "LVGU_058.fastq.gz"
mv "LG_GP_barcode03_merged_porechop-abi.fastq.gz" "LVGU_059.fastq.gz"
mv "merged_NB01.fastq.gz" "LVGU_001.fastq.gz"
mv "merged_NB03-NB02.fastq.gz" "LVGU_004.fastq.gz"
mv "merged_NB04-NB03.fastq.gz" "LVGU_005.fastq.gz"
mv "merged_NB06-NB05.fastq.gz" "LVGU_010.fastq.gz"
mv "merged_NB08-NB06.fastq.gz" "LVGU_012.fastq.gz"
mv "merged_NB09-NB11.fastq.gz" "LVGU_013.fastq.gz"
mv "merged_NB10-NB12.fastq.gz" "LVGU_014.fastq.gz"
mv "merged_NB12-NB04.fastq.gz" "LVGU_006.fastq.gz"
mv "merged_NB13.fastq.gz" "LVGU_018.fastq.gz"
mv "merged_NB14.fastq.gz" "LVGU_020.fastq.gz"
mv "merged_NB15.fastq.gz" "LVGU_023.fastq.gz"
mv "merged_NB16.fastq.gz" "LVGU_024.fastq.gz"
mv "single_NB02.fastq.gz" "LVGU_002.fastq.gz"
mv "single_NB05.fastq.gz" "LVGU_009.fastq.gz"
mv "single_NB07.fastq.gz" "LVGU_011.fastq.gz"
mv "single_NB11.fastq.gz" "LVGU_015.fastq.gz"
mv "single_NB17.fastq.gz" "LVGU_025.fastq.gz"

echo "Renaming complete!"
echo "Files renamed:"
ls -1 LVGU_*.fastq.gz | sort -V
