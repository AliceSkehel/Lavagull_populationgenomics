# ONT vs Illumina SNP Comparison Table
library(seqinr)

# Set working directory
setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Load sequences
illumina_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ont_seqs <- read.fasta("ont_all_angsd.fasta", as.string = TRUE)

# Clean up sample names - remove "LVGU_" prefix and standardize numbers
names(illumina_seqs) <- gsub("LVGU_", "", names(illumina_seqs))
names(ont_seqs) <- gsub("LVGU_", "", names(ont_seqs))

# Convert ONT sample names to remove leading zeros
names(ont_seqs) <- as.character(as.numeric(names(ont_seqs)))

# Also standardize Illumina names to ensure they're just numbers
names(illumina_seqs) <- as.character(as.numeric(names(illumina_seqs)))

# Find individuals present in both datasets
common_individuals <- intersect(names(illumina_seqs), names(ont_seqs))

cat("Found", length(common_individuals), "individuals in both datasets:\n")
cat(paste(common_individuals, collapse = ", "), "\n\n")

if(length(common_individuals) == 0) {
  cat("No common individuals found! Check sample naming.\n")
  cat("Illumina samples:", paste(names(illumina_seqs)[1:5], collapse = ", "), "...\n")
  cat("ONT samples:", paste(names(ont_seqs)[1:5], collapse = ", "), "...\n")
} else {
  
  # Function to compare two sequences and count differences
  compare_sequences <- function(seq1, seq2) {
    chars1 <- strsplit(seq1, "")[[1]]
    chars2 <- strsplit(seq2, "")[[1]]
    
    # Make sure sequences are same length
    min_length <- min(length(chars1), length(chars2))
    chars1 <- chars1[1:min_length]
    chars2 <- chars2[1:min_length]
    
    # Count positions where both have data (not N)
    valid_positions <- (chars1 != "N" & chars1 != "n") & (chars2 != "N" & chars2 != "n")
    
    if(sum(valid_positions) == 0) {
      return(list(differences = NA, total_compared = 0, percent_diff = NA))
    }
    
    # Count differences at valid positions
    differences <- sum(chars1[valid_positions] != chars2[valid_positions])
    total_compared <- sum(valid_positions)
    percent_diff <- (differences / total_compared) * 100
    
    return(list(differences = differences, total_compared = total_compared, percent_diff = percent_diff))
  }
  
  # Create comparison table
  comparison_results <- data.frame(
    Individual = character(0),
    SNP_Differences = numeric(0),
    Sites_Compared = numeric(0),
    Percent_Different = numeric(0),
    Illumina_Coverage = numeric(0),
    ONT_Coverage = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Compare each individual
  for(individual in common_individuals) {
    illumina_seq <- illumina_seqs[[individual]]
    ont_seq <- ont_seqs[[individual]]
    
    # Calculate coverage for each technology
    illumina_chars <- strsplit(illumina_seq, "")[[1]]
    ont_chars <- strsplit(ont_seq, "")[[1]]
    
    illumina_coverage <- sum(illumina_chars != "N" & illumina_chars != "n") / length(illumina_chars)
    ont_coverage <- sum(ont_chars != "N" & ont_chars != "n") / length(ont_chars)
    
    # Compare sequences
    comparison <- compare_sequences(illumina_seq, ont_seq)
    
    # Add to results
    comparison_results <- rbind(comparison_results, data.frame(
      Individual = individual,
      SNP_Differences = comparison$differences,
      Sites_Compared = comparison$total_compared,
      Percent_Different = round(comparison$percent_diff, 3),
      Illumina_Coverage = round(illumina_coverage * 100, 1),
      ONT_Coverage = round(ont_coverage * 100, 1),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by number of differences
  comparison_results <- comparison_results[order(comparison_results$SNP_Differences), ]
  
  # Print results
  cat("=== ONT vs ILLUMINA COMPARISON RESULTS ===\n\n")
  print(comparison_results, row.names = FALSE)
  
  # Summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  valid_comparisons <- !is.na(comparison_results$SNP_Differences)
  
  if(sum(valid_comparisons) > 0) {
    cat("Mean SNP differences:", round(mean(comparison_results$SNP_Differences, na.rm = TRUE), 2), "\n")
    cat("Median SNP differences:", median(comparison_results$SNP_Differences, na.rm = TRUE), "\n")
    cat("Range of differences:", range(comparison_results$SNP_Differences, na.rm = TRUE), "\n")
    cat("Mean percent different:", round(mean(comparison_results$Percent_Different, na.rm = TRUE), 3), "%\n")
    
    # Identify best and worst concordance
    best_match <- comparison_results[which.min(comparison_results$SNP_Differences), ]
    worst_match <- comparison_results[which.max(comparison_results$SNP_Differences), ]
    
    cat("\nBest concordance:", best_match$Individual, "with", best_match$SNP_Differences, "differences\n")
    cat("Worst concordance:", worst_match$Individual, "with", worst_match$SNP_Differences, "differences\n")
    
    # Coverage comparison
    cat("\nMean Illumina coverage:", round(mean(comparison_results$Illumina_Coverage), 1), "%\n")
    cat("Mean ONT coverage:", round(mean(comparison_results$ONT_Coverage), 1), "%\n")
  }
  
}
 
