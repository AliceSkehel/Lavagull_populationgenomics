# Analyze Evolutionary Relationships Between Haplotypes (with Ancient DNA)
library(ape)
library(seqinr)

setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Load temporal haplotype data (includes ancient samples)
haplotypes <- read.csv("analysis_outputs/temporal_haplotypes_clustered.csv")

# Load modern and ancient sequences
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)
seqs <- c(modern_seqs, ancient_seqs)

# Filter by coverage
coverage <- sapply(seqs, function(x) {
  chars <- strsplit(x, "")[[1]]
  sum(chars != "N" & chars != "n") / length(chars)
})
seqs <- seqs[coverage >= 0.3]

# Get one representative per haplotype
hap_reps <- sapply(unique(haplotypes$Haplotype), function(h) {
  samples <- haplotypes$Sample[haplotypes$Haplotype == h]
  samples[1]
})

# Calculate distances between haplotype representatives
rep_seqs <- seqs[hap_reps]
rep_matrix <- do.call(rbind, lapply(rep_seqs, function(x) strsplit(x, "")[[1]]))
rownames(rep_matrix) <- names(rep_seqs)
rep_matrix[rep_matrix == "n"] <- "N"
rep_dna <- as.DNAbin(rep_matrix)
hap_distances <- dist.dna(rep_dna, model = "raw", pairwise.deletion = TRUE)

# Convert distances to matrix for easier reading
dist_matrix <- as.matrix(hap_distances)
hap_labels <- haplotypes$Haplotype[match(rownames(dist_matrix), haplotypes$Sample)]
rownames(dist_matrix) <- hap_labels
colnames(dist_matrix) <- hap_labels

cat("=== HAPLOTYPE DISTANCE MATRIX ===\n")
print(round(dist_matrix, 4))

# Calculate summary statistics
all_distances <- dist_matrix[upper.tri(dist_matrix)]
cat("Minimum distance between haplotypes:", round(min(all_distances), 4), "\n")
cat("Maximum distance between haplotypes:", round(max(all_distances), 4), "\n")
cat("Mean distance between haplotypes:", round(mean(all_distances), 4), "\n")

# Simple tree to visualize relationships
hap_tree <- nj(hap_distances)
cat("Tree created with", length(hap_tree$tip.label), "haplotypes\n")

# Save the distance matrix
write.csv(dist_matrix, "analysis_outputs/haplotype_distances.csv")