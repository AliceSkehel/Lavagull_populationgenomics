#!/bin/bash

# Input BEAGLE file
BEAGLE=~/Sequencing_Combined/angsd_variant_calling_all35/lava_gulls_all35_merged.beagle.gz

# Output directory
OUTDIR=~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis
mkdir -p "$OUTDIR"

echo "Running NGSadmix for K=1 to K=10, 6 replicates each with different seeds..."
echo "Start time: $(date)"

for K in {1..10}; do
    echo "Running K=$K..."
    for REP in {1..6}; do
        # Set different seeds for each replicate
        SEED=$((REP * 250000))
        echo "  K=$K, Replicate $REP, Seed=$SEED..."
        NGSadmix -likes "$BEAGLE" -seed $SEED -K "$K" -o "$OUTDIR/lava_gulls_k${K}_rep${REP}" -P 8 -minMaf 0.05
    done
    echo "All replicates for K=$K complete!"
done

echo "All K values (1-10) complete!"
echo "End time: $(date)"
