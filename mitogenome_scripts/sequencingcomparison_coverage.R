# mtDNA Coverage Analysis - Minimal Version
library(seqinr)
library(ggplot2)
library(reshape2)
setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Read sequences and calculate coverage
illumina_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE, forceDNAtolower = FALSE)
ont_seqs <- read.fasta("ont_all_angsd.fasta", as.string = TRUE, forceDNAtolower = FALSE)

seq_length <- nchar(illumina_seqs[[1]])
cat(paste("Analyzing", seq_length, "bp mtDNA genome\n"))

# Calculate coverage per position (% samples with non-N bases)
get_coverage <- function(seqs) {
  coverage_matrix <- matrix(0, nrow = length(seqs), ncol = seq_length)
  for(i in 1:length(seqs)) {
    chars <- strsplit(seqs[[i]], "")[[1]]
    coverage_matrix[i, ] <- ifelse(chars == "N", 0, 1)
  }
  return(colMeans(coverage_matrix) * 100)
}

illumina_cov <- get_coverage(illumina_seqs)
ont_cov <- get_coverage(ont_seqs)

# Create coverage data
coverage_data <- data.frame(
  Position = 1:seq_length,
  Illumina = illumina_cov,
  ONT = ont_cov,
  Average = (illumina_cov + ont_cov) / 2
)

# Plot coverage
coverage_long <- melt(coverage_data[,1:3], id.vars = "Position", 
                      variable.name = "Platform", value.name = "Coverage")

p1 <- ggplot(coverage_long, aes(x = Position, y = Coverage, color = Platform)) +
  geom_line(alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE) +
  scale_color_manual(values = c("Illumina" = "blue", "ONT" = "red")) +
  labs(title = "mtDNA Positional Coverage", x = "Position (bp)", 
       y = "Coverage (% samples)", color = "Platform") +
  theme_minimal() + ylim(0, 100)

# Plot difference
p2 <- ggplot(coverage_data, aes(x = Position, y = Illumina - ONT)) +
  geom_line(color = "purple", alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Coverage Difference (Illumina - ONT)", x = "Position (bp)", 
       y = "Difference (%)") + theme_minimal()

# Display plots
print(p1)
print(p2)

# Summary statistics
cat(paste("Mean coverage - Illumina:", round(mean(illumina_cov), 1), 
          "%, ONT:", round(mean(ont_cov), 1), "%\n"))

# Find worst regions (consecutive positions <50% coverage)
find_bad_regions <- function(coverage, threshold = 50) {
  bad_pos <- which(coverage < threshold)
  if(length(bad_pos) == 0) return(NULL)
  
  regions <- c()
  start <- bad_pos[1]
  for(i in 2:length(bad_pos)) {
    if(bad_pos[i] != bad_pos[i-1] + 1) {
      if(bad_pos[i-1] - start >= 9) {  # Regions ≥10bp
        regions <- rbind(regions, c(start, bad_pos[i-1], bad_pos[i-1] - start + 1))
      }
      start <- bad_pos[i]
    }
  }
  if(bad_pos[length(bad_pos)] - start >= 9) {
    regions <- rbind(regions, c(start, bad_pos[length(bad_pos)], 
                                bad_pos[length(bad_pos)] - start + 1))
  }
  if(!is.null(regions)) colnames(regions) <- c("Start", "End", "Length")
  return(regions)
}

ill_bad <- find_bad_regions(illumina_cov)
ont_bad <- find_bad_regions(ont_cov)

cat("\nProblem regions (≥10bp with <50% coverage):\n")
if(!is.null(ill_bad)) { cat("Illumina:\n"); print(ill_bad) }
if(!is.null(ont_bad)) { cat("ONT:\n"); print(ont_bad) }

# Show 10 worst positions
worst_pos <- order(coverage_data$Average)[1:10]
cat("\n10 worst positions:\n")
print(coverage_data[worst_pos, c("Position", "Average", "Illumina", "ONT")])