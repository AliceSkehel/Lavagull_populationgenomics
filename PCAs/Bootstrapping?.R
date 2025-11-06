
pacman::p_load(boot)


# Only the numeric SNP columns
beagle <- read.table("lava_gulls_all35_merged.beagle.gz", header=TRUE)
geno_mat <- as.matrix(beagle[,-1])  # remove first column (positions)
rownames(geno_mat) <- beagle[,1]    # optional: SNP IDs
