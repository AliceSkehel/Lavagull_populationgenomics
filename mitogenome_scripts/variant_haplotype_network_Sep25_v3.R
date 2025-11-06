# SELF-CONTAINED mtDNA Dendrogram Script
library(seqinr)
library(dendextend)
library(paletteer)

# Set working directory (adjust path as needed)
setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Load modern and ancient sequences
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)

# Combine sequences
seqs <- c(modern_seqs, ancient_seqs)

# Clean up sample names - remove "LVGU_" prefix
names(seqs) <- gsub("LVGU_", "", names(seqs))

# Calculate coverage and filter
coverage <- sapply(seqs, function(x) {
  chars <- strsplit(x, "")[[1]]
  sum(chars != "N" & chars != "n") / length(chars)
})
seqs <- seqs[coverage >= 0.6]  # Your 60% coverage threshold

# Create sequence matrix
seq_matrix <- do.call(rbind, lapply(seqs, function(x) strsplit(x, "")[[1]]))
rownames(seq_matrix) <- names(seqs)

# Fix missing data notation
seq_matrix[seq_matrix == "n"] <- "N"

# Find variable positions
variable_positions <- c()
for(i in 1:ncol(seq_matrix)) {
  col <- seq_matrix[, i]
  valid_bases <- col[col != "N"]
  if(length(unique(valid_bases)) > 1 & length(valid_bases) > 0) {
    variable_positions <- c(variable_positions, i)
  }
}
snp_matrix <- seq_matrix[, variable_positions]

# Calculate pairwise distances
calculate_pairwise_distance <- function(snp_matrix) {
  n_samples <- nrow(snp_matrix)
  dist_matrix <- matrix(0, n_samples, n_samples)
  rownames(dist_matrix) <- rownames(snp_matrix)
  colnames(dist_matrix) <- rownames(snp_matrix)
  
  for(i in 1:n_samples) {
    for(j in 1:n_samples) {
      if(i != j) {
        seq1 <- snp_matrix[i, ]
        seq2 <- snp_matrix[j, ]
        
        valid_pos <- (seq1 != "N") & (seq2 != "N")
        
        if(sum(valid_pos) > 0) {
          differences <- sum(seq1[valid_pos] != seq2[valid_pos])
          total_compared <- sum(valid_pos)
          dist_matrix[i, j] <- differences / total_compared
        } else {
          dist_matrix[i, j] <- 1
        }
      }
    }
  }
  
  return(as.dist(dist_matrix))
}

# Calculate distances
snp_distances <- calculate_pairwise_distance(snp_matrix)

# Create hierarchical clustering
distances <- as.dist(snp_distances)
hc <- hclust(distances, method = "average")
dend <- as.dendrogram(hc)

# Get the fishualize palette colors
fish_colors <- paletteer_d("fishualize::Acanthurus_coeruleus", 4)

# Create color vector - all black by default
branch_colors <- rep("black", 13)

# Set clusters 10-13 with fishualize palette colors
branch_colors[10] <- fish_colors[1]  # Cluster 10
branch_colors[11] <- fish_colors[2]  # Cluster 11
branch_colors[12] <- fish_colors[3]  # Cluster 12
branch_colors[13] <- fish_colors[4]  # Cluster 13

# Set margins and plot
par(mar=c(9,1,1,1))

options(axes = TRUE)

# Create the beautiful dendrogram
dend %>%
  set("labels_col", value = branch_colors, k=13) %>%
  set("branches_k_color", value = branch_colors, k=13) %>%
  plot(main = "mtDNA Phylogeny - Hierarchical Clustering", 
       horiz=FALSE, axes=FALSE)

text(x = 0, y = c(0, 0.2, 0.4, 0.6), labels = c("0", "0.2", "0.4", "0.6"), 
     pos = 2, xpd = TRUE)
axis(2)
mtext("Genetic Distance (proportion of variable sites)", side = 2, line = 2.5, cex = 1.1)

# Add legend
legend("topright", 
       legend = c("Cluster 10", "Cluster 11", "Cluster 12", "Cluster 13", "Other clusters"),
       col = c(fish_colors[1], fish_colors[2], fish_colors[3], fish_colors[4], "black"),
       lty = 1, lwd = 2, cex = 0.8, bty = "n", title = "Highlighted Clusters")

cat("Self-contained dendrogram complete!\n")
cat("Found", length(variable_positions), "variable positions\n")
cat("Analyzing", nrow(snp_matrix), "samples\n")

