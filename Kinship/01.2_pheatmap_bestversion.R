# Load required packages
pacman::p_load(tidyverse, pheatmap, pvclust, reshape2)

# Set working directory
setwd("~/Sequencing_Combined/angsd_variant_calling_all35/")

# Load site/metadata
Kinship_Summary_Data <- read.csv("Kinship_Summary_Data.csv", header=TRUE)


kinship <- read.table("lava_gulls_kinship_clean.txt", header=TRUE)

# Read sample mapping
samples <- read.table("sample_mapping.txt", header=FALSE, sep="\t")
colnames(samples) <- c("Ind", "Sample")

# View dimensions
cat("Total pairwise comparisons:", nrow(kinship), "\n")
cat("Number of samples:", nrow(samples), "\n")

cat("\n=== KING coefficient summary ===\n")
summary(kinship$KING)

# Filter out low coverage ancient samples and different sequencing technique
samples_to_exclude <- c("LVGU_473", "LVGU_476", "LVGU_60")

# Get their Ind numbers
exclude_ind <- samples$Ind[samples$Sample %in% samples_to_exclude]
exclude_ind_num <- as.numeric(gsub("Ind", "", exclude_ind))

cat("Excluding samples:", samples_to_exclude, "\n")
cat("Ind numbers:", exclude_ind_num, "\n")

# Filter kinship data (remove pairs involving excluded samples)
kinship_filtered <- kinship[!(kinship$ida %in% exclude_ind_num | 
                                kinship$idb %in% exclude_ind_num), ]

cat("\nOriginal pairs:", nrow(kinship), "\n")
cat("Filtered pairs:", nrow(kinship_filtered), "\n")
cat("Samples remaining:", length(unique(c(kinship_filtered$ida, kinship_filtered$idb))), "\n")

# Now recreate the matrix with only filtered samples
remaining_samples <- sort(unique(c(kinship_filtered$ida, kinship_filtered$idb)))
n_remaining <- length(remaining_samples)

rab_matrix_filtered <- matrix(NA, nrow=n_remaining, ncol=n_remaining)

# Map old indices to new matrix positions
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

# Get sample names for remaining samples
remaining_sample_names <- samples$Sample[match(remaining_samples, 
                                               as.numeric(gsub("Ind", "", samples$Ind)))]
rownames(rab_matrix_filtered) <- remaining_sample_names
colnames(rab_matrix_filtered) <- remaining_sample_names


kinship <- read.table("lava_gulls_kinship_clean.txt", header=TRUE)

# Standardize ID column to uppercase for matching
Kinship_Summary_Data$ID <- toupper(Kinship_Summary_Data$ID)

# Create annotation dataframe with Island, Site, Sex, and Batch
annotation_df <- data.frame(
  row.names = remaining_sample_names,
  Island = Kinship_Summary_Data$`ISLAND.or.Country`[match(remaining_sample_names, Kinship_Summary_Data$ID)],
  Site = Kinship_Summary_Data$Locality[match(remaining_sample_names, Kinship_Summary_Data$ID)],
  Sex = Kinship_Summary_Data$Sex[match(remaining_sample_names, Kinship_Summary_Data$ID)],
  Batch = Kinship_Summary_Data$Transit.Batches[match(remaining_sample_names, Kinship_Summary_Data$ID)]
) %>% 
  view()


# Hierarchical clustering with bootstrap using pvclust
set.seed(123)
pv_clustering <- pvclust(as.matrix(rab_matrix_filtered), 
                         method.hclust = "average",
                         method.dist = "euclidean", 
                         nboot = 10000)  # This will take a while!

# Define annotation colors for the color bars
annotation_colors <- list(
  Island = c("Genovesa" = "#E69F00",       # Orange
             "San Cristóbal" = "#56B4E9",   # Blue  
             "Santa Cruz" = "#009E73",
             "Isabela" = "purple"),     # Green
  Sex = c("Male" = "#0072B2", 
          "Female" = "#D55E00", 
          "Unknown" = "gray80")
)

# Create heatmap with pvclust clustering
pheatmap(rab_matrix_filtered,
         cluster_rows = pv_clustering$hclust,
         cluster_cols = pv_clustering$hclust,
         color = color_palette,
         main = "",
         scale = "none",
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 50,
         treeheight_col = 50)


# Create heatmap WITH annotation bars
pheatmap(rab_matrix_filtered,
         cluster_rows = pv_clustering$hclust,
         cluster_cols = pv_clustering$hclust,
         annotation_row = annotation_df,      # Add color bars on rows
         annotation_col = annotation_df,      # Add color bars on columns
         annotation_colors = annotation_colors,  # Define colors
         color = color_palette,
         main = "",
         scale = "none",
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 50,
         treeheight_col = 50)


# Create annotation with Island as fallback for missing Sites
annotation_df <- data.frame(
  row.names = remaining_sample_names,
  Island = Kinship_Summary_Data$ISLAND.or.Country[match(remaining_sample_names, Kinship_Summary_Data$ID)],
  Site = Kinship_Summary_Data$Locality[match(remaining_sample_names, Kinship_Summary_Data$ID)]
)

# For samples with missing Site, use Island instead
annotation_df$Site[is.na(annotation_df$Site)] <- annotation_df$Island[is.na(annotation_df$Site)]

# Create simplified annotation
annotation_df_simple <- data.frame(
  row.names = remaining_sample_names,
  Site = annotation_df$Site
)

# Colors with Santiago and Fernandina added
annotation_colors <- list(
  Site = c(
    # Ancient sample islands (add new colors)
    "Santiago" = "#E31A1C",      # Red
    "Fernandina" = "#FF7F00",    # Orange-red
    
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
# Create color palette (they used black-white, you can customize)
color_palette <- colorRampPalette(c("white", "grey30", "black"))(100)


# Create heatmap with only Site annotation
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
         breaks = seq(0.15, 1, length.out = 101))  # Emphasize your actual range

