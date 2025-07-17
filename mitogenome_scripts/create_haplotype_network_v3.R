#!/usr/bin/env Rscript

# Set up user library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
}
.libPaths(c(user_lib, .libPaths()))

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Load required libraries
suppressMessages({
  library(ape)
  library(pegas)
  library(phangorn)
})

# Set working directory
setwd("/home/askehel/mitogenome/haplotype_analysis")

cat("=== Lava Gull Haplotype Network Analysis ===\n")

# Check if file exists
if (!file.exists("aligned_sequences.fasta")) {
  stop("Error: aligned_sequences.fasta not found!")
}

# Read aligned sequences
tryCatch({
  sequences <- read.dna("aligned_sequences.fasta", format = "fasta")
  cat("Successfully read", nrow(sequences), "sequences\n")
  cat("Sequence length:", ncol(sequences), "bp\n\n")
}, error = function(e) {
  stop("Error reading FASTA file: ", e$message)
})

# Filter out reference sequence (keep only Lava Gull sequences)
seq_names <- rownames(sequences)
lava_gull_indices <- grep("LG|Lava|barcode", seq_names, ignore.case = TRUE)

if (length(lava_gull_indices) == 0) {
  cat("Warning: No Lava Gull sequences found. Using all sequences.\n")
  lava_sequences <- sequences
} else {
  lava_sequences <- sequences[lava_gull_indices, ]
  cat("Filtered to", nrow(lava_sequences), "Lava Gull sequences\n")
}

# Convert to haplotypes
cat("Identifying haplotypes...\n")
hap <- haplotype(lava_sequences)
cat("Number of unique haplotypes:", length(hap), "\n\n")

# Print haplotype assignments
hap_assignments <- attr(hap, "index")
for (i in 1:length(hap)) {
  individuals <- which(hap_assignments == i)
  cat("Haplotype", i, ":", rownames(lava_sequences)[individuals], "\n")
}

# Create haplotype network
cat("\nCreating haplotype network...\n")
net <- haploNet(hap)

# Create network plot
png("lava_gull_haplotype_network.png", width = 1000, height = 800, res = 150)
plot(net, size = attr(net, "freq"), scale.ratio = 2, 
     cex = 1.2, pie = table(attr(lava_sequences, "dimnames")[[1]]),
     bg = "lightblue", labels = TRUE)
title("Lava Gull Mitochondrial Haplotype Network", cex.main = 1.5)
dev.off()

# Calculate distance matrix
cat("Calculating genetic distances...\n")
dist_matrix <- dist.dna(lava_sequences, model = "raw")

# Create summary
cat("\n=== Genetic Distance Summary ===\n")
cat("Mean pairwise distance:", round(mean(dist_matrix), 6), "\n")
cat("Maximum distance:", round(max(dist_matrix), 6), "\n")
cat("Minimum distance:", round(min(dist_matrix[dist_matrix > 0]), 6), "\n")

# Save results
write.csv(as.matrix(dist_matrix), "pairwise_distances.csv")

# Create detailed summary file
sink("haplotype_analysis_results.txt")
cat("=== Lava Gull Haplotype Analysis Results ===\n")
cat("Analysis date:", date(), "\n")
cat("Total sequences analyzed:", nrow(lava_sequences), "\n")
cat("Sequence length:", ncol(lava_sequences), "bp\n")
cat("Unique haplotypes found:", length(hap), "\n")
cat("Mean pairwise distance:", mean(dist_matrix), "\n")
cat("Maximum pairwise distance:", max(dist_matrix), "\n")
cat("Minimum pairwise distance:", min(dist_matrix[dist_matrix > 0]), "\n\n")

cat("Haplotype assignments:\n")
for (i in 1:length(hap)) {
  individuals <- which(attr(hap, "index") == i)
  cat("Haplotype", i, ":", rownames(lava_sequences)[individuals], "\n")
}
sink()

cat("\n=== Analysis Complete! ===\n")
cat("Files created:\n")
cat("- lava_gull_haplotype_network.png (network visualization)\n")
cat("- pairwise_distances.csv (genetic distances)\n")
cat("- haplotype_analysis_results.txt (detailed summary)\n")
