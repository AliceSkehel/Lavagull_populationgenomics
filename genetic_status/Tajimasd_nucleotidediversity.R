setwd("~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis")

# Read merged data with full path
data <- read.table("genome_wide_thetas_combined.txt", header=FALSE, skip=1)

# Manually assign correct column names
colnames(data) <- c("indexInfo", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

# Now check
head(data)

# Calculate weighted averages
total_sites <- sum(data$nSites)
genome_tW <- sum(data$tW * data$nSites) / total_sites
genome_tP <- sum(data$tP * data$nSites) / total_sites
genome_Tajima <- sum(data$Tajima * data$nSites) / total_sites

# Print results
cat("\n==============================================\n")
cat("GENOME-WIDE DIVERSITY STATISTICS\n")
cat("(33 modern lava gull samples)\n")
cat("==============================================\n\n")
cat("Total sites analyzed:", formatC(total_sites, format="d", big.mark=","), "\n\n")
cat("Watterson's theta (tW):", formatC(genome_tW, format="f", digits=2, big.mark=","), "\n")
cat("Nucleotide diversity (π):", formatC(genome_tP, format="f", digits=2, big.mark=","), "\n")
cat("Tajima's D:", round(genome_Tajima, 3), "\n\n")

cat("Interpretation:\n")
cat("- Tajima's D =", round(genome_Tajima, 2), "(negative = excess rare variants)\n")
cat("- Suggests: population expansion or purifying selection\n")
cat("- Low π indicates reduced genetic diversity\n")
cat("==============================================\n\n")
