library(seqinr)
library(ape)

setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")
dir.create("analysis_outputs", showWarnings = FALSE)

# Read modern and ancient sequences
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)

# Keep original names and combine
seqs <- c(modern_seqs, ancient_seqs)

# Filter by coverage (lower threshold for ancient DNA)
coverage <- sapply(seqs, function(x) {
  chars <- strsplit(x, "")[[1]]
  sum(chars != "N" & chars != "n") / length(chars)
})
good_samples <- coverage >= 0.3
seqs <- seqs[good_samples]
cat("Kept", length(seqs), "samples with >=30% coverage\n")

# Convert to DNA format and calculate distances
seq_matrix <- do.call(rbind, lapply(seqs, function(x) strsplit(x, "")[[1]]))
rownames(seq_matrix) <- names(seqs)
seq_matrix[seq_matrix == "n"] <- "N"
dna <- as.DNAbin(seq_matrix)
distances <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE)

# Cluster by similarity (0.1% divergence threshold)
clusters <- cutree(hclust(as.dist(distances)), h = 0.001)

# Create results with time period info
results <- data.frame(
  Sample = names(seqs),
  Haplotype = paste0("H", clusters),
  Time_Period = ifelse(names(seqs) %in% c("LVGU_473", "LVGU_476"), "Ancient", "Modern")
)

# Save and print
cat("Found", length(unique(clusters)), "haplotypes using similarity clustering\n")
print(sort(table(results$Haplotype), decreasing = TRUE))

# Print assignments by time period
for(period in c("Ancient", "Modern")) {
  cat("\n", period, " samples:\n", sep = "")
  period_data <- results[results$Time_Period == period, ]
  for(hap in sort(unique(period_data$Haplotype))) {
    samples <- period_data$Sample[period_data$Haplotype == hap]
    cat(hap, "(n=", length(samples), "): ", paste(samples, collapse = ", "), "\n", sep = "")
  }
}

# Plot dendrogram
plot(hclust(as.dist(distances)), main = "Temporal Clustering of Sequences")
abline(h = 0.001, col = "red")

