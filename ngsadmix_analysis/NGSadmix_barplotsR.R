library(tidyverse)
library(RColorBrewer)

# Set working directory
setwd("~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis/")

# Read sample names
samples <- read.table("sample_names.txt", header = FALSE)$V1

# Read best replicates based on our summary
k2 <- read.table("lava_gulls_k2_rep4.qopt")
k3 <- read.table("lava_gulls_k3_rep2.qopt")
k4 <- read.table("lava_gulls_k4_rep5.qopt")
k5 <- read.table("lava_gulls_k5_rep5.qopt")

# Add sample names to rows
rownames(k2) <- samples
rownames(k3) <- samples
rownames(k4) <- samples
rownames(k5) <- samples

# Quick check
print(dim(k2))
print(head(samples))

# Read the metadata (e.g., islands)
metadata <- read.csv("~/Sequencing_Combined/angsd_variant_calling_all35/Kinship_Summary_Data.csv")

# Create a vector of island names for each sample
island_metadata <- metadata %>%
  filter(ISLAND.or.Country != "") %>%  # Make sure there is an island listed
  dplyr::select(ID, ISLAND.or.Country)

# Match samples with their respective islands
sample_island_mapping <- match(samples, island_metadata$ID)
sample_islands <- island_metadata$ISLAND.or.Country[sample_island_mapping]

# Set color palette for K values (i.e., the lineage colors)
cluster_colors_k2 <- RColorBrewer::brewer.pal(2, "Pastel1")  # 2 colors for K=2
cluster_colors_k3 <- RColorBrewer::brewer.pal(3, "Pastel1")  # 3 colors for K=3
cluster_colors_k4 <- RColorBrewer::brewer.pal(4, "Pastel1")  # 4 colors for K=4
cluster_colors_k5 <- RColorBrewer::brewer.pal(5, "Pastel1")  # 5 colors for K=5

# Reorder samples based on island to group island samples together
ordered_samples <- samples[order(sample_islands)]

# Reorder K data to match the ordered samples (so the bar plots reflect island clusters)
k2_ordered <- k2[match(ordered_samples, rownames(k2)),]
k3_ordered <- k3[match(ordered_samples, rownames(k3)),]
k4_ordered <- k4[match(ordered_samples, rownames(k4)),]
k5_ordered <- k5[match(ordered_samples, rownames(k5)),]

# Set up plotting layout (4 plots in a single column)
par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))

# Function to plot each K value with the island grouping and color by K (lineage)
plot_island_grouped <- function(k_data, cluster_colors, ordered_samples, main_title) {
  barplot(t(as.matrix(k_data)), 
          col = rep(cluster_colors, length.out = ncol(k_data)),  # Color by K lineage
          border = NA,
          space = 0,
          xlab = ,
          ylab = "Ancestry",
          main = main_title,
          names.arg = ordered_samples,
          las = 2,
          cex.names = 0.7)
}

# Plot for K=2 (with lineage colors)
plot_island_grouped(k2_ordered, cluster_colors_k2, ordered_samples, "K=2")

# Plot for K=3 (with lineage colors)
plot_island_grouped(k3_ordered, cluster_colors_k3, ordered_samples, "K=3")

# Plot for K=4 (with lineage colors)
plot_island_grouped(k4_ordered, cluster_colors_k4, ordered_samples, "K=4")

# Plot for K=5 (with lineage colors)
plot_island_grouped(k5_ordered, cluster_colors_k5, ordered_samples, "K=5")

