#!/usr/bin/env Rscript

# Haplotype Network Analysis for Lava Gull mtDNA
# Install required packages if needed
if (!require("ape")) install.packages("ape")
if (!require("pegas")) install.packages("pegas")
if (!require("phangorn")) install.packages("phangorn")

library(ape)
library(pegas)
library(phangorn)

# Set working directory
setwd("/home/askehel/mitogenome/haplotype_analysis")

# Read aligned sequences
sequences <- read.dna("aligned_sequences.fasta", format = "fasta")

# Print basic information
cat("=== Lava Gull Haplotype Network Analysis ===\n")
cat("Number of sequences:", nrow(sequences), "\n")
cat("Sequence length:", ncol(sequences), "bp\n\n")

# Convert to haplotypes
hap <- haplotype(sequences)
cat("Number of unique haplotypes:", length(hap), "\n")

# Create haplotype network
net <- haploNet(hap)

# Save the network plot
png("lava_gull_haplotype_network.png", width = 800, height = 600)
plot(net, size = attr(net, "freq"), scale.ratio = 2, 
     cex = 0.8, pie = table(attr(sequences, "dimnames")[[1]]))
title("Lava Gull Mitochondrial Haplotype Network")
dev.off()

# Create a distance matrix
dist_matrix <- dist.dna(sequences, model = "raw")
cat("\n=== Pairwise Distances ===\n")
print(as.matrix(dist_matrix))

# Save detailed results
write.csv(as.matrix(dist_matrix), "pairwise_distances.csv")

# Create a simple summary
sink("haplotype_summary.txt")
cat("=== Lava Gull Haplotype Analysis Summary ===\n")
cat("Total sequences analyzed:", nrow(sequences), "\n")
cat("Unique haplotypes found:", length(hap), "\n")
cat("Average pairwise distance:", mean(dist_matrix), "\n")
cat("Maximum pairwise distance:", max(dist_matrix), "\n")
cat("Minimum pairwise distance:", min(dist_matrix), "\n")
sink()

cat("Analysis complete! Check these files:\n")
cat("- lava_gull_haplotype_network.png (network visualization)\n")
cat("- pairwise_distances.csv (genetic distances)\n")
cat("- haplotype_summary.txt (summary statistics)\n")

