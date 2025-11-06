#!/bin/bash

# Input BEAGLE file
BEAGLE=~/Sequencing_Combined/angsd_variant_calling_all35/lava_gulls_all35_merged.beagle.gz

# Output directory (make sure it exists)
OUTDIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
mkdir -p "$OUTDIR"

echo "Running NGSadmix for K=2, 6 replicates..."
echo "Start time: $(date)"

K=2

for REP in {1..6}; do
    echo "  K=$K, Replicate $REP..."
    NGSadmix -likes "$BEAGLE" -K "$K" -o "$OUTDIR/lava_gulls_k${K}_rep${REP}" -P 8 -minMaf 0.05
done

echo "All replicates for K=$K complete!"
echo "End time: $(date)"

