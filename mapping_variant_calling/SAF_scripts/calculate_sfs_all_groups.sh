#!/bin/bash
WORK_DIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis

cd $WORK_DIR
echo "Step 2: Calculating SFS from merged file..."
realSFS lava_gulls_all_merged.saf.idx > lava_gulls_all_merged.sfs
