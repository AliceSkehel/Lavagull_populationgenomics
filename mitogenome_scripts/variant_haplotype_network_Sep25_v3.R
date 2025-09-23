# Use paletteer Antique palette for better colors
library(dendextend)
library(paletteer)

distances <- as.dist(snp_distances)
hc <- hclust(distances, method = "average")
dend <- as.dendrogram(hc)

# Get the fishualize Acanthurus_coeruleus palette colors
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

dend %>%
  set("labels_col", value = branch_colors, k=13) %>%
  set("branches_k_color", value = branch_colors, k=13) %>%
  plot(main = "mtDNA Phylogeny - Hierarchical Clustering", 
       horiz=FALSE, axes=FALSE)

axis(2)
mtext("Genetic Distance (proportion of variable sites)", side = 2, line = 2.5, cex = 1.1)

legend("topright", 
       legend = c("Cluster 10", "Cluster 11", "Cluster 12", "Cluster 13", "Other clusters"),
       col = c(fish_colors[1], fish_colors[2], fish_colors[3], fish_colors[4], "black"),
       lty = 1, lwd = 2, cex = 0.8, bty = "n", title = "Highlighted Clusters")

cat("Cluster 13 = light blue, Cluster 10 = light green, rest = black!\n")