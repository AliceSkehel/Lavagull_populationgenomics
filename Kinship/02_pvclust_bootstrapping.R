# Load required packages
pacman::p_load(tidyverse, pheatmap, pvclust, reshape2)

# Set working directory
setwd("~/Sequencing_Combined/angsd_variant_calling_all35/")

# Load data
Kinship_Summary_Data <- read.csv("Kinship_Summary_Data.csv", header=TRUE)
kinship <- read.table("lava_gulls_kinship_clean.txt", header=TRUE)
samples <- read.table("sample_mapping.txt", header=FALSE, sep="\t")
colnames(samples) <- c("Ind", "Sample")

# View dimensions
summary(kinship$KING)

# Filter out low coverage ancient samples and different sequencing technique
samples_to_exclude <- c("LVGU_473", "LVGU_476", "LVGU_60")

exclude_ind <- samples$Ind[samples$Sample %in% samples_to_exclude]
exclude_ind_num <- as.numeric(gsub("Ind", "", exclude_ind))

# Filter kinship data
kinship_filtered <- kinship[!(kinship$ida %in% exclude_ind_num | 
                                kinship$idb %in% exclude_ind_num), ]

# Create matrix with filtered samples (32 modern, high-coverage samples)
remaining_samples <- sort(unique(c(kinship_filtered$ida, kinship_filtered$idb)))
n_remaining <- length(remaining_samples)
rab_matrix_filtered <- matrix(NA, nrow=n_remaining, ncol=n_remaining)
idx_map <- setNames(1:n_remaining, remaining_samples)

for(i in 1:nrow(kinship_filtered)) {
  ida <- kinship_filtered$ida[i]
  idb <- kinship_filtered$idb[i]
  rab_val <- kinship_filtered$rab[i]
  
  new_i <- idx_map[as.character(ida)]
  new_j <- idx_map[as.character(idb)]
  
  rab_matrix_filtered[new_i, new_j] <- rab_val
  rab_matrix_filtered[new_j, new_i] <- rab_val
}

diag(rab_matrix_filtered) <- 1

# Add sample names
remaining_sample_names <- samples$Sample[match(remaining_samples, 
                                               as.numeric(gsub("Ind", "", samples$Ind)))]
rownames(rab_matrix_filtered) <- remaining_sample_names
colnames(rab_matrix_filtered) <- remaining_sample_names

# Prepare annotations
Kinship_Summary_Data$ID <- toupper(Kinship_Summary_Data$ID)

annotation_df <- data.frame(
  row.names = remaining_sample_names,
  Island = Kinship_Summary_Data$ISLAND.or.Country[match(remaining_sample_names, Kinship_Summary_Data$ID)],
  Site = Kinship_Summary_Data$Locality[match(remaining_sample_names, Kinship_Summary_Data$ID)]
)

# Create simplified annotation with Site only
annotation_df_simple <- data.frame(
  row.names = remaining_sample_names,
  Site = annotation_df$Site
)

# Bootstrap clustering (THIS WILL TAKE TIME!)
set.seed(123)
pv_clustering <- pvclust(as.matrix(rab_matrix_filtered), 
                         method.hclust = "average",
                         method.dist = "euclidean", 
                         nboot = 10000)

# Site colors
annotation_colors <- list(
  Site = c(
    # Genovesa sites (orange)
    "Darwins Bay" = "#E69F00",
    
    # San Cristóbal sites (shades of blue)
    "Loberia" = "#0072B2",
    "Playa de Oro" = "#4A9BD9",
    "Playa de los Marinos" = "#5FA8D3",
    "Rosa Blanca" = "#7FB3E0",
    "Playa Ochoa" = "#99C4E8",
    "Bahía Sardinas" = "#B4D7F0",
    
    # Santa Cruz sites (shades of green)
    "Charles Darwin Station Pier" = "#009E73",
    "Tortuga Bay" = "#4DB89A",
    "Playa de los Alemanas" = "#66C4A8",
    "Fishermans pier" = "#80D1B8",
    
    # Isabela sites (shades of purple)
    "Malecon" = "#9370DB",
    "Muelle" = "#A98DE3",
    "El Faro/Surf Spot" = "#BFA9E9",
    "Playa del Amor" = "#D5C6EF"
  )
)

# Color palette
color_palette <- colorRampPalette(c("white", "grey30", "black"))(100)

# Create final heatmap
pheatmap(rab_matrix_filtered,
         cluster_rows = pv_clustering$hclust,
         cluster_cols = pv_clustering$hclust,
         annotation_row = annotation_df_simple,
         annotation_col = annotation_df_simple,
         annotation_colors = annotation_colors,
         color = color_palette,
         main = "",
         scale = "none",
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 60,
         treeheight_col = 60,
         legend = FALSE,
         breaks = seq(0.15, 1, length.out = 101))

#### Create R0-R1 scatterplot ####

par(mfrow=c(1,2))  # 1 row, 2 columns

plot(kinship_filtered$R1, kinship_filtered$R0,
     xlab="R1", ylab="R0",
     pch=19, col=rgb(0,0,0,0.3),
     cex=0.5)

# Identify parent-offspring pairs
po_candidates <- kinship_filtered[kinship_filtered$R0 < 0.1 & 
                                    kinship_filtered$R1 > 0.9, ]

if(nrow(po_candidates) > 0) {
  points(po_candidates$R1, po_candidates$R0, col="red", pch=19, cex=1)
}

# Create KING-R1 scatterplot
plot(kinship_filtered$R1, kinship_filtered$KING,
     xlab="R1", ylab="KING",
     pch=19, col=rgb(0,0,0,0.3),
     cex=0.5)

if(nrow(po_candidates) > 0) {
  points(po_candidates$R1, po_candidates$KING, col="red", pch=19, cex=1)
}

threshold_line <- (min(first_degree_KING) + max(non_first_degree_KING)) / 2
abline(h=threshold_line, lty=2, col="blue")


# Print inferred parent-offspring pairs
for(i in 1:nrow(po_candidates)) {
  sample_a <- samples$Sample[samples$Ind == paste0("Ind", po_candidates$ida[i])]
  sample_b <- samples$Sample[samples$Ind == paste0("Ind", po_candidates$idb[i])]
  cat(sample_a, "vs", sample_b,
      "| R0=", round(po_candidates$R0[i], 3),
      "| R1=", round(po_candidates$R1[i], 3),
      "| KING=", round(po_candidates$KING[i], 3), "\n")
}


# Add Island to your annotation
annotation_df_simple$Island <- annotation_df$Island

# Test within-island vs between-island relatedness
observed_within_island <- mean(rab_matrix_filtered[
  outer(annotation_df_simple$Island, annotation_df_simple$Island, "==") & 
    lower.tri(rab_matrix_filtered)
], na.rm=TRUE)

observed_between_island <- mean(rab_matrix_filtered[
  outer(annotation_df_simple$Island, annotation_df_simple$Island, "!=") & 
    lower.tri(rab_matrix_filtered)
], na.rm=TRUE)

# Permutation test (island level)
n_permutations <- 10000
perm_diffs_island <- numeric(n_permutations)

for(i in 1:n_permutations) {
  shuffled_islands <- sample(annotation_df_simple$Island)
  
  within_perm <- mean(rab_matrix_filtered[
    outer(shuffled_islands, shuffled_islands, "==") & 
      lower.tri(rab_matrix_filtered)
  ], na.rm=TRUE)
  
  between_perm <- mean(rab_matrix_filtered[
    outer(shuffled_islands, shuffled_islands, "!=") & 
      lower.tri(rab_matrix_filtered)
  ], na.rm=TRUE)
  
  perm_diffs_island[i] <- within_perm - between_perm
}

p_value_island <- mean(perm_diffs_island >= (observed_within_island - observed_between_island))
cat("Permutation test p-value:", p_value_island, "\n")

if(p_value_island < 0.05) {
  cat("SIGNIFICANT: Individuals from same island are more related\n")
} else {
  cat("NOT SIGNIFICANT: No island-level structure\n")
}

# Plot
hist(perm_diffs_island, breaks=50, main="Permutation Test - Island Level", 
     xlab="Difference in mean relatedness (within - between island)")
abline(v=observed_within_island - observed_between_island, col="red", lwd=2)

# Show island distribution
cat("\nIsland distribution:\n")
print(table(annotation_df_simple$Island))

