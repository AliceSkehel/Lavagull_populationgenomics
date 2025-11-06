
pacman::p_load(ggplot2, ggforce, dplyr, tidyverse, patchwork)
setwd("~/Sequencing_Combined/angsd_variant_calling_all35/")

site_colors <- c(
  "Santiago" = "#E31A1C",
  "Fernandina" = "#FF7F00",
  "Darwins Bay" = "#E69F00",
  "Loberia" = "#0072B2",
  "Playa de Oro" = "#4A9BD9",
  "Playa de los Marinos" = "#5FA8D3",
  "Rosa Blanca" = "#7FB3E0",
  "Playa Ochoa" = "#99C4E8",
  "Bahía Sardinas" = "#B4D7F0",
  "Charles Darwin Station Pier" = "#009E73",
  "Tortuga Bay" = "#4DB89A",
  "Playa de los Alemanas" = "#66C4A8",
  "Fishermans pier" = "#80D1B8",
  "Malecon" = "#9370DB",
  "Muelle" = "#A98DE3",
  "El Faro/Surf Spot" = "#BFA9E9",
  "Playa del Amor" = "#D5C6EF"
)

# Load covariance matrix and PCA
cov_matrix1 <- as.matrix(read.table("pcangsd/lava_gulls_filtered29_e2.cov"))
pca1 <- eigen(cov_matrix1)
scores1 <- pca$vectors[, 1:2]

# Load metadata
metadata <- read.csv("Kinship_Summary_Data.csv", header=TRUE)

# Filter to keep only the 27 samples (remove IDs 2, 8, 9, 13, 23, 25, 33, 34)
# Removed: LVGU_3, LVGU_10, LvGu_11, LvGu_15, LVGU_39, LVGU_42, LVGU_473, LVGU_476
removed_samples1 <- c("LVGU_3", "LVGU_10", "LvGu_11", "LvGu_15", "LVGU_39", "LVGU_42")
metadata_filtered1 <- metadata %>% 
  filter(!ID %in% removed_samples) %>%
  filter(ID != "LVGU_60")  # Also remove LVGU_60 (reference genome)

# Assign sample IDs to PCA scores
rownames(scores1) <- metadata_filtered1$ID

# Combine into a dataframe
pca_df2 <- data.frame(
  SampleID = rownames(scores1),
  PC1 = scores[,1],
  PC2 = scores[,2]) %>% 
  left_join(metadata_filtered1, by = c("SampleID" = "ID"))


# Plot
p3_all <- ggplot(pca_df2, aes(x = PC1, y = PC2, color = Locality, shape = ISLAND.or.Country)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Locality), type = "t", level = 0.95, linetype = 2, size = 1) +
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pca$values[1] / sum(pca$values) * 100, 1), "%)"),
       y = paste0("PC2 (", round(pca$values[2] / sum(pca$values) * 100, 1), "%)"),
       title = NULL, color = NULL, shape = NULL) +
  theme(text = element_text(size = 14), legend.title = element_blank())


p5 + p3_all + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')



# Prepare the data for plotting
variance_df <- data.frame(
  PC = paste0("PC", 1:length(pca$values)), 
  Variance_explained = (pca$values / sum(pca$values)) * 100
)

# Sort by variance in descending order and ensure PC is a factor
variance_df <- variance_df %>% 
  arrange(desc(Variance_explained))
variance_df$PC <- factor(variance_df$PC, levels = variance_df$PC)

# Plot the eigenvalues (variance explained)
scree_plot <- ggplot(variance_df, aes(x = PC, y = Variance_explained)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_line(aes(group = 1), color = "#9370DB", size = 1) +
  geom_point(color = "#9370DB", size = 2) +
  labs(
    title = NULL,
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(scree_plot)
