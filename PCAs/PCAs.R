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
cov_matrix <- as.matrix(read.table("lava_gulls_all35_PCA.cov"))
pca <- eigen(cov_matrix)
scores <- pca$vectors[, 1:2]

# Load metadata and remove the rogue sample
metadata <- read.csv("Kinship_Summary_Data.csv", header=TRUE)
metadata_filtered <- metadata[metadata$ID != "LVGU_53", ]

# Assign sample IDs to PCA scores
rownames(scores) <- metadata_filtered$ID

# Combine into a dataframe
pca_df <- data.frame(
  SampleID = rownames(scores),
  PC1 = scores[,1],
  PC2 = scores[,2]
) %>% 
  left_join(metadata_filtered, by = c("SampleID" = "ID"))

p5 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Locality, shape = ISLAND.or.Country)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Locality), type = "t", level = 0.95, linetype = 2, size = 1) +
  scale_color_manual(values = site_colors) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(pca$values[1] / sum(pca$values) * 100, 1), "%)"),
    y = paste0("PC2 (", round(pca$values[2] / sum(pca$values) * 100, 1), "%)"),
    title = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme(
    text = element_text(size = 14),
    legend.title = element_blank()
  )
print(p5)
