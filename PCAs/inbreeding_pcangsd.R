# Inbreeding pcangsd

pacman::p_load(ggplot2,dplyr,tidyverse)

setwd("~/Sequencing_Combined/angsd_variant_calling_all35/")

# Read inbreeding coefficients
inbreeding <- read.table("pcangsd/lava_gulls_all35_inbreeding.inbreed.samples", header=FALSE)
colnames(inbreeding) <- "F"

# Add sample IDs
metadata <- read.csv("Kinship_Summary_Data.csv", header=TRUE)

# Filter metadata to match the 35 samples in your analysis
metadata_filtered <- metadata %>% 
  filter(ID != "LVGU_53")      # Remove LVGU_53

inbreeding$Sample <- metadata_filtered$ID

# Merge with metadata
inbreeding <- left_join(inbreeding, metadata_filtered, by=c("Sample"="ID"))

# View the data
head(inbreeding)
summary(inbreeding$F)

island_colors <- c(
  "Santiago" = "#E57373",
  "Fernandina" = "#FF7F00",
  "Genovesa" = "#E69F00",
  "San Cristóbal" = "#5FA8D3",
  "Santa Cruz" = "#4DB89A",
  "Isabela" = "#9370DB"
)

# Plot by island with new colors
p1 <- ggplot(inbreeding, aes(x=ISLAND.or.Country, y=F, fill=ISLAND.or.Country)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.2, size=3, alpha=0.6) +
  scale_fill_manual(values=island_colors) +
  theme_bw() +
  labs(y="Inbreeding Coefficient (F)", 
       x="Island",
       title=NULL) +
  theme(text = element_text(size=14),
        legend.position = "none")

# Plot by individual with island colors
p2 <- ggplot(inbreeding, aes(x=reorder(Sample, F), y=F, fill=ISLAND.or.Country)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=island_colors) +
  coord_flip() +
  theme_bw() +
  labs(x="Sample", 
       y="Inbreeding Coefficient (F)",
       fill="Island",
       title=NULL) +
  theme(text = element_text(size=12))

print(p1)
print(p2)

# Check the unique values in your island column
unique(inbreeding$ISLAND.or.Country)

# If they don't match exactly, you might need to recode them
# For example, if it's "San Cristobal" instead of "San Cristóbal":
inbreeding <- inbreeding %>%
  mutate(Island_fixed = case_when(
    ISLAND.or.Country == "San Cristobal" ~ "San Cristóbal",
    TRUE ~ ISLAND.or.Country
  ))

# Then use the fixed column
p2 <- ggplot(inbreeding, aes(x=reorder(Sample, F), y=F, fill=Island_fixed)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=island_colors) +
  coord_flip() +
  theme_bw() +
  labs(x="Sample", 
       y="Inbreeding Coefficient (F)",
       fill="Island",
       title=NULL) +
  theme(text = element_text(size=12))

print(p2)

