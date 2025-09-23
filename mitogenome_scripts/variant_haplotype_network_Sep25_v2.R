# Beautiful mtDNA Phylogenetic Tree Visualization
library(ggtree)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ape)

# Assuming you have nj_tree and sample_info from previous analysis

# Enhanced sample information with better aesthetics
sample_info_enhanced <- sample_info %>%
  mutate(
    Time_Period_Label = case_when(
      Time_Period == "Ancient" ~ "Historical (1905-1932)",
      Time_Period == "Modern" ~ "Contemporary",
      TRUE ~ Time_Period
    ),
    Haplotype_Group = case_when(
      Haplotype == "V2" ~ "Major Modern Lineage",
      Haplotype %in% c("V12", "V13") ~ "Historical Lineages", 
      TRUE ~ "Minor Modern Lineages"
    ),
    Point_Size = case_when(
      Time_Period == "Ancient" ~ 4,
      Haplotype == "V2" ~ 3,
      TRUE ~ 2.5
    )
  )

# Define a beautiful color palette
time_colors <- c("Historical (1905-1932)" = "#E31A1C", 
                 "Contemporary" = "#1F78B4")

haplotype_colors <- c("Historical Lineages" = "#E31A1C",
                      "Major Modern Lineage" = "#1F78B4", 
                      "Minor Modern Lineages" = "#33A02C")

# Check the actual tree depth
cat("Maximum tree depth:", max(node.depth.edgelength(nj_tree)), "\n")
cat("Maximum pairwise distance:", max(snp_distances), "\n")

# Fix negative edge lengths in the tree
nj_tree$edge.length[nj_tree$edge.length < 0] <- 0

# Create a haplotype-focused plot
p3 <- ggtree(nj_tree, layout = "rectangular", linewidth = 0.8, color = "grey40") +
  geom_tippoint(aes(color = sample_info_enhanced$Haplotype_Group[match(label, sample_info_enhanced$Sample)],
                    size = sample_info_enhanced$Point_Size[match(label, sample_info_enhanced$Sample)]),
                alpha = 0.8) +
  geom_tiplab(aes(label = paste0(label, " (", 
                                 sample_info_enhanced$Haplotype[match(label, sample_info_enhanced$Sample)], ")")), 
              size = 2.8, 
              hjust = -0.1, 
              fontface = "italic",
              color = "black") +
  scale_color_manual(values = haplotype_colors, name = "Lineage Group") +
  scale_size_identity() +
  theme_tree2(panel.grid.major = element_line(color = "grey95", size = 0.2),
              panel.grid.minor = element_blank()) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  labs(title = "",
       subtitle = "") +
  xlim(0, max(node.depth.edgelength(nj_tree)) * 1.6)

print(p3)


