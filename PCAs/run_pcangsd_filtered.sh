#!/bin/bash

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pcangsd_env

# Run pcangsd
cd ~/Sequencing_Combined/angsd_variant_calling_all35/pcangsd/

pcangsd.py -beagle ../lava_gulls_filtered27.beagle.gz \
           -o lava_gulls_filtered27 \
           -threads 4

# Convert covariance matrix to text
python -c "import numpy as np; cov_matrix = np.load('lava_gulls_filtered27.cov.npy'); np.savetxt('lava_gulls_filtered27.cov', cov_matrix, fmt='%0.10f')"

echo "PCAngsd analysis complete!"
