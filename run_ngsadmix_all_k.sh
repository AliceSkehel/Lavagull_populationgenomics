#!/bin/bash

BEAGLE=~/Sequencing_Combined/angsd_variant_calling_all35/lava_gulls_all35_merged.beagle.gz
OUTDIR=~/Lavagull_populationgenomics/ngsadmix_analysis

cd $OUTDIR

echo "Running NGSadmix for K=2 to K=6..."
echo "Start time: $(date)"

# K=2
echo "Running K=2..."
NGSadmix -likes $BEAGLE -K 2 -o lava_gulls_k2 -P 8 -minMaf 0.05

# K=3
echo "Running K=3..."
NGSadmix -likes $BEAGLE -K 3 -o lava_gulls_k3 -P 8 -minMaf 0.05

# K=4
echo "Running K=4..."
NGSadmix -likes $BEAGLE -K 4 -o lava_gulls_k4 -P 8 -minMaf 0.05

# K=5
echo "Running K=5..."
NGSadmix -likes $BEAGLE -K 5 -o lava_gulls_k5 -P 8 -minMaf 0.05

# K=6
echo "Running K=6..."
NGSadmix -likes $BEAGLE -K 6 -o lava_gulls_k6 -P 8 -minMaf 0.05

echo "All K values complete! End time: $(date)"
