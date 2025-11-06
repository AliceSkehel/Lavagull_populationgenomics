# ASkehel
# 28th Oct
# Kinship attempt #1 - following the killer whale kinship paper using NgsRelate
# Have 35 samples and x2 ancient
# used angsd script saved in /home/askehel/Lavagull_populationgenomics/mapping_variant_calling/angsd_"g1-4c"

pacman::p_load(reshape2, pheatmap, pvclust, tidyverse)

setwd("~/Sequencing_Combined/angsd_variant_calling_all35/")

kinship <- read.table("lava_gulls_kinship_clean.txt", header=TRUE)

# Read sample mapping
samples <- read.table("sample_mapping.txt", header=FALSE, sep="\t")
colnames(samples) <- c("Ind", "Sample")

# View dimensions
cat("Total pairwise comparisons:", nrow(kinship), "\n")
cat("Number of samples:", nrow(samples), "\n")

cat("\n=== KING coefficient summary ===\n")
summary(kinship$KING)

cat("\n=== Close relatives (KING > 0.177 = 1st degree) ===\n")
close_rels <- kinship[kinship$KING > 0.177, ]
close_rels_sorted <- close_rels[order(-close_rels$KING), c("ida","idb","KING","R0","R1")]
print(close_rels_sorted)

# Map to sample names
close_rels_sorted$Sample_A <- samples$Sample[match(close_rels_sorted$ida, 
                                                   gsub("Ind","",samples$Ind))]
close_rels_sorted$Sample_B <- samples$Sample[match(close_rels_sorted$idb, 
                                                   gsub("Ind","",samples$Ind))]
print(close_rels_sorted[,c("Sample_A","Sample_B","KING","R0","R1")])

cat("\n=== Ancient samples comparison ===\n")
ancient_pair <- kinship[(kinship$ida==33 & kinship$idb==34) | 
                          (kinship$ida==34 & kinship$idb==33), ]
print(ancient_pair[, c("ida","idb","KING","R0","R1")])
cat("LVGU_473 vs LVGU_476 are 1st degree relatives!\n")

# Create distance matrix from rab (relatedness coefficient)
# Using 1-rab as distance metric (like the paper)
kinship_wide <- acast(kinship, ida ~ idb, value.var="rab")

# Convert to distance matrix
dist_matrix <- as.dist(1 - kinship_wide)

# Hierarchical clustering with UPGMA
cat("\n=== Performing hierarchical clustering ===\n")
hc <- hclust(dist_matrix, method="average")

# Plot dendrogram directly to screen
plot(hc, main="Lava Gull Kinship Dendrogram (UPGMA)", 
     xlab="Individual", ylab="Distance (1-relatedness)")

# Create heatmap directly to screen
pheatmap(kinship_wide, 
         main="Lava Gull Pairwise Relatedness (rab)",
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         clustering_method="average")


# Create a symmetric matrix from pairwise kinship data
# First, create empty matrix
n_samples <- 35
rab_matrix <- matrix(NA, nrow=n_samples, ncol=n_samples)
rownames(rab_matrix) <- 0:(n_samples-1)
colnames(rab_matrix) <- 0:(n_samples-1)

# Fill in the matrix (both upper and lower triangles)
for(i in 1:nrow(kinship)) {
  ida <- kinship$ida[i]
  idb <- kinship$idb[i]
  rab_val <- kinship$rab[i]
  
  # Fill both triangles
  rab_matrix[ida+1, idb+1] <- rab_val
  rab_matrix[idb+1, ida+1] <- rab_val
}

# Diagonal should be 1 (individual with itself)
diag(rab_matrix) <- 1

# Check for NAs
cat("Number of NAs in matrix:", sum(is.na(rab_matrix)), "\n")

# Convert to distance matrix
dist_matrix <- as.dist(1 - rab_matrix)

# Now clustering should work
hc <- hclust(dist_matrix, method="average")

# Plot dendrogram
plot(hc, main="Lava Gull Kinship Dendrogram (UPGMA)", 
     xlab="Individual", ylab="Distance (1-relatedness)")

# Heatmap
pheatmap(rab_matrix, 
         main="Lava Gull Pairwise Relatedness (rab)",
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         clustering_method="average")

sample_names <- samples$Sample[match(0:(n_samples-1), 
                                     as.numeric(gsub("Ind", "", samples$Ind)))]

# Apply sample names to matrix
rownames(rab_matrix) <- sample_names
colnames(rab_matrix) <- sample_names

# Check it worked
cat("Sample names applied:\n")
print(head(rownames(rab_matrix)))

# Now when you plot, it will show actual sample names
# Convert to distance matrix
dist_matrix <- as.dist(1 - rab_matrix)

# Clustering
hc <- hclust(dist_matrix, method="average")

# Plot dendrogram with sample names
plot(hc, main="Lava Gull Kinship Dendrogram (UPGMA)", 
     xlab="Sample", ylab="Distance (1-relatedness)",
     cex=0.6)  # Make text smaller if needed

# Heatmap with sample names
pheatmap(rab_matrix, 
         main="Lava Gull Pairwise Relatedness (rab)",
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         clustering_method="average",
         fontsize_row=8,  # Adjust font size
         fontsize_col=8)

