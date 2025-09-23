# Variant-Based Temporal Haplotype Analysis
library(seqinr)

setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Load modern and ancient sequences
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)

# Combine sequences (keep original names)
seqs <- c(modern_seqs, ancient_seqs)

# Calculate coverage and filter
coverage <- sapply(seqs, function(x) {
  chars <- strsplit(x, "")[[1]]
  sum(chars != "N" & chars != "n") / length(chars)
})
seqs <- seqs[coverage >= 0.6]

seq_matrix <- do.call(rbind, lapply(seqs, function(x) strsplit(x, "")[[1]]))
rownames(seq_matrix) <- names(seqs)

# Fix missing data notation
cat("Before fixing: unique characters in matrix:", paste(unique(as.vector(seq_matrix))[1:10], collapse = ", "), "\n")
seq_matrix[seq_matrix == "n"] <- "N"
cat("After fixing: unique characters in matrix:", paste(unique(as.vector(seq_matrix))[1:10], collapse = ", "), "\n")

# Find variable positions
variable_positions <- c()
for(i in 1:ncol(seq_matrix)) {
  col <- seq_matrix[, i]
  valid_bases <- col[col != "N"]
  if(length(unique(valid_bases)) > 1 & length(valid_bases) > 0) {
    variable_positions <- c(variable_positions, i)
  }
}
snp_matrix <- seq_matrix[, variable_positions]

cat("Found", length(variable_positions), "variable positions\n")

# Compare sequences function (handles missing data)
compare_sequences_fixed <- function(seq1, seq2, min_overlap = 50) {  # Lower for ancient DNA
  
  valid_positions <- (seq1 != "N") & (seq2 != "N")
  
  if(sum(valid_positions) < min_overlap) {
    return(NA)
  }
  
  identical_at_valid <- all(seq1[valid_positions] == seq2[valid_positions])
  return(identical_at_valid)
}

# Group sequences function
group_sequences_fixed <- function(snp_matrix, min_overlap = 50) {
  
  n_samples <- nrow(snp_matrix)
  sample_names <- rownames(snp_matrix)
  
  groups <- as.list(1:n_samples)
  names(groups) <- sample_names
  
  cat("Comparing", n_samples, "samples...\n")
  
  for(i in 1:(n_samples-1)) {
    for(j in (i+1):n_samples) {
      
      seq1 <- snp_matrix[i, ]
      seq2 <- snp_matrix[j, ]
      sample1 <- sample_names[i]
      sample2 <- sample_names[j]
      
      identical <- compare_sequences_fixed(seq1, seq2, min_overlap)
      
      if(!is.na(identical) && identical) {
        group1 <- groups[[sample1]]
        group2 <- groups[[sample2]]
        
        merged_group <- unique(c(group1, group2))
        
        for(sample in sample_names[merged_group]) {
          groups[[sample]] <- merged_group
        }
        
        cat("Grouped", sample1, "and", sample2, "\n")
      }
    }
  }
  
  # Extract unique groups
  unique_groups <- unique(groups)
  group_assignments <- rep(NA, n_samples)
  names(group_assignments) <- sample_names
  
  for(i in seq_along(unique_groups)) {
    group_members <- sample_names[unique_groups[[i]]]
    group_assignments[group_members] <- i
  }
  
  return(group_assignments)
}

# Run the variant-based grouping
cat("\nRunning variant-based temporal grouping...\n")
variant_groups <- group_sequences_fixed(snp_matrix, min_overlap = 50)

# Create results with time period info
variant_results <- data.frame(
  Sample = names(variant_groups),
  Haplotype = paste0("V", variant_groups),
  Time_Period = ifelse(names(variant_groups) %in% c("LVGU_473", "LVGU_476"), "Ancient", "Modern"),
  stringsAsFactors = FALSE
)

# Show results
hap_freq_variant <- table(variant_results$Haplotype)
cat("\nVariant-based haplotype frequencies:\n")
print(sort(hap_freq_variant, decreasing = TRUE))

# Print assignments by time period
for(period in c("Ancient", "Modern")) {
  cat("\n", period, " samples:\n", sep = "")
  period_data <- variant_results[variant_results$Time_Period == period, ]
  for(hap in sort(unique(period_data$Haplotype))) {
    samples <- period_data$Sample[period_data$Haplotype == hap]
    cat(hap, "(n=", length(samples), "): ", paste(samples, collapse = ", "), "\n", sep = "")
  }
}

# Analyze temporal patterns
cat("\n=== TEMPORAL LINEAGE ANALYSIS ===\n")
ancient_haps <- unique(variant_results$Haplotype[variant_results$Time_Period == "Ancient"])
modern_haps <- unique(variant_results$Haplotype[variant_results$Time_Period == "Modern"])
shared_haps <- intersect(ancient_haps, modern_haps)

cat("Ancient haplotypes:", paste(ancient_haps, collapse = ", "), "\n")
cat("Modern haplotypes:", paste(modern_haps, collapse = ", "), "\n")
cat("Shared haplotypes:", ifelse(length(shared_haps) > 0, paste(shared_haps, collapse = ", "), "None"), "\n")

# Calculate lineage survival
if(length(ancient_haps) > 0) {
  survival_rate <- length(shared_haps) / length(ancient_haps) * 100
  cat("Lineage survival rate:", round(survival_rate, 1), "%\n")
}

# Test different overlap requirements
cat("\n=== TESTING DIFFERENT OVERLAP REQUIREMENTS ===\n")
for(min_overlap in c(25, 50, 100, 200)) {
  test_groups <- group_sequences_fixed(snp_matrix, min_overlap = min_overlap)
  n_groups <- length(unique(test_groups))
  cat("Min overlap", min_overlap, "bp: ", n_groups, "haplotype groups\n")
}

# Save results
write.csv(variant_results, "analysis_outputs/variant_temporal_haplotypes.csv", row.names = FALSE)

cat("\nVariant-based temporal analysis complete!\n")
cat("Files saved:\n")
cat("- variant_temporal_haplotypes.csv\n")