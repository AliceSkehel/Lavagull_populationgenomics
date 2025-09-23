# MtDNA kinship relatedness
# Alice
# 2025-07-17
# using aligned_v4.fasta

library(ape)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggtree)

# Read your aligned sequences
sequences <- read.dna("/home/askehel/Sequencing_Combined/mtdna_only/consensus/aligned_v4.fasta", format = "fasta")

# Basic info about your dataset
cat("Number of sequences:", nrow(sequences), "\n")
cat("Sequence length:", ncol(sequences), "\n")
cat("Sample names:", rownames(sequences), "\n")

# Calculate genetic distances
dist_matrix <- dist.dna(sequences, model = "TN93")

# Convert to full matrix
kinship_matrix <- as.matrix(dist_matrix)

# Convert distance to relatedness (inverse)
relatedness_matrix <- 1 - kinship_matrix

# Create phylogenetic tree
nj_tree <- nj(dist_matrix)

# Create combined tree + heatmap plot
pheatmap(relatedness_matrix, 
         color = colorRampPalette(c("white", "lightblue", "steelblue", "darkblue"))(100),
         main = "Maternal Kinship: Tree + Heatmap",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         clustering_distance_rows = dist_matrix,
         clustering_distance_cols = dist_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE)


pheatmap(relatedness_matrix, 
         color = colorRampPalette(c("white", "lightblue", "steelblue", "darkblue"))(100),
         main = "Maternal Kinship: Tree + Heatmap",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         clustering_distance_rows = dist_matrix,
         clustering_distance_cols = dist_matrix,
         clustering_method = "ward.D2",  # Better clustering method
         cluster_rows = TRUE,
         cluster_cols = TRUE)

