# Comprehensive ONT vs Illumina Comparison with Multiple References
library(seqinr)

# Set working directory
setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Load sequences
illumina_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ont_seqs <- read.fasta("ont_all_angsd.fasta", as.string = TRUE)

# Load the yellow-legged gull reference genome
ref_genome <- read.fasta("mtdna_reference_only.fna", as.string = TRUE)
yllg_reference_seq <- ref_genome[[1]]

# Clean up sample names
names(illumina_seqs) <- gsub("LVGU_", "", names(illumina_seqs))
names(ont_seqs) <- gsub("LVGU_", "", names(ont_seqs))
names(ont_seqs) <- as.character(as.numeric(names(ont_seqs)))
names(illumina_seqs) <- as.character(as.numeric(names(illumina_seqs)))

# Get individual 60 as reference
individual_60_ref <- NULL
if("60" %in% names(illumina_seqs)) {
  individual_60_ref <- illumina_seqs[["60"]]
  cat("Using Illumina individual 60 as second reference\n")
} else if("60" %in% names(ont_seqs)) {
  individual_60_ref <- ont_seqs[["60"]]
  cat("Using ONT individual 60 as second reference\n")
} else {
  cat("Individual 60 not found - will only use yellow-legged gull reference\n")
}

# Find common individuals
common_individuals <- intersect(names(illumina_seqs), names(ont_seqs))
common_individuals <- common_individuals[common_individuals != "60"]  # Remove reference individual

cat("Found", length(common_individuals), "individuals for comparison:\n")
cat(paste(common_individuals, collapse = ", "), "\n\n")

# Function to compare sequence against reference
compare_to_reference <- function(seq, ref_seq) {
  chars <- strsplit(seq, "")[[1]]
  ref_chars <- strsplit(ref_seq, "")[[1]]
  
  min_length <- min(length(chars), length(ref_chars))
  chars <- chars[1:min_length]
  ref_chars <- ref_chars[1:min_length]
  
  valid_positions <- (chars != "N" & chars != "n") & (ref_chars != "N" & ref_chars != "n")
  
  if(sum(valid_positions) == 0) {
    return(list(differences = NA, total_compared = 0, percent_diff = NA))
  }
  
  differences <- sum(chars[valid_positions] != ref_chars[valid_positions])
  total_compared <- sum(valid_positions)
  percent_diff <- (differences / total_compared) * 100
  
  return(list(differences = differences, total_compared = total_compared, percent_diff = percent_diff))
}

# Function to compare two sequences directly
compare_sequences <- function(seq1, seq2) {
  chars1 <- strsplit(seq1, "")[[1]]
  chars2 <- strsplit(seq2, "")[[1]]
  
  min_length <- min(length(chars1), length(chars2))
  chars1 <- chars1[1:min_length]
  chars2 <- chars2[1:min_length]
  
  valid_positions <- (chars1 != "N" & chars1 != "n") & (chars2 != "N" & chars2 != "n")
  
  if(sum(valid_positions) == 0) {
    return(list(differences = NA, total_compared = 0, percent_diff = NA))
  }
  
  differences <- sum(chars1[valid_positions] != chars2[valid_positions])
  total_compared <- sum(valid_positions)
  percent_diff <- (differences / total_compared) * 100
  
  return(list(differences = differences, total_compared = total_compared, percent_diff = percent_diff))
}

# TABLE 1: ONT vs Illumina Direct Comparison
cat("=== TABLE 1: ONT vs ILLUMINA DIRECT COMPARISON ===\n\n")

ont_vs_illumina_results <- data.frame(
  Individual = character(0),
  SNP_Differences = numeric(0),
  Sites_Compared = numeric(0),
  Percent_Different = numeric(0),
  Illumina_Coverage = numeric(0),
  ONT_Coverage = numeric(0),
  stringsAsFactors = FALSE
)

for(individual in common_individuals) {
  if(individual %in% names(illumina_seqs) && individual %in% names(ont_seqs)) {
    illumina_seq <- illumina_seqs[[individual]]
    ont_seq <- ont_seqs[[individual]]
    
    # Calculate coverage
    illumina_chars <- strsplit(illumina_seq, "")[[1]]
    ont_chars <- strsplit(ont_seq, "")[[1]]
    illumina_coverage <- sum(illumina_chars != "N" & illumina_chars != "n") / length(illumina_chars)
    ont_coverage <- sum(ont_chars != "N" & ont_chars != "n") / length(ont_chars)
    
    # Skip if coverage too low
    if(illumina_coverage < 0.6 | ont_coverage < 0.6) {
      cat("Skipping individual", individual, "- Low coverage\n")
      next
    }
    
    # Compare sequences
    comparison <- compare_sequences(illumina_seq, ont_seq)
    
    ont_vs_illumina_results <- rbind(ont_vs_illumina_results, data.frame(
      Individual = individual,
      SNP_Differences = comparison$differences,
      Sites_Compared = comparison$total_compared,
      Percent_Different = round(comparison$percent_diff, 3),
      Illumina_Coverage = round(illumina_coverage * 100, 1),
      ONT_Coverage = round(ont_coverage * 100, 1),
      stringsAsFactors = FALSE
    ))
  }
}

ont_vs_illumina_results <- ont_vs_illumina_results[order(ont_vs_illumina_results$SNP_Differences), ]
print(ont_vs_illumina_results, row.names = FALSE)

# TABLE 2: Comparison against Yellow-legged Gull Reference
cat("\n\n=== TABLE 2: COMPARISON AGAINST YELLOW-LEGGED GULL REFERENCE ===\n\n")

yllg_ref_results <- data.frame(
  Individual = character(0),
  Illumina_vs_YLLG_Diffs = numeric(0),
  Illumina_vs_YLLG_Sites = numeric(0),
  Illumina_vs_YLLG_Percent = numeric(0),
  ONT_vs_YLLG_Diffs = numeric(0),
  ONT_vs_YLLG_Sites = numeric(0),
  ONT_vs_YLLG_Percent = numeric(0),
  Illumina_Coverage = numeric(0),
  ONT_Coverage = numeric(0),
  stringsAsFactors = FALSE
)

for(individual in common_individuals) {
  if(individual %in% names(illumina_seqs) && individual %in% names(ont_seqs)) {
    illumina_seq <- illumina_seqs[[individual]]
    ont_seq <- ont_seqs[[individual]]
    
    # Calculate coverage
    illumina_chars <- strsplit(illumina_seq, "")[[1]]
    ont_chars <- strsplit(ont_seq, "")[[1]]
    illumina_coverage <- sum(illumina_chars != "N" & illumina_chars != "n") / length(illumina_chars)
    ont_coverage <- sum(ont_chars != "N" & ont_chars != "n") / length(ont_chars)
    
    if(illumina_coverage < 0.6 | ont_coverage < 0.6) next
    
    # Compare to YLLG reference
    illumina_vs_yllg <- compare_to_reference(illumina_seq, yllg_reference_seq)
    ont_vs_yllg <- compare_to_reference(ont_seq, yllg_reference_seq)
    
    yllg_ref_results <- rbind(yllg_ref_results, data.frame(
      Individual = individual,
      Illumina_vs_YLLG_Diffs = illumina_vs_yllg$differences,
      Illumina_vs_YLLG_Sites = illumina_vs_yllg$total_compared,
      Illumina_vs_YLLG_Percent = round(illumina_vs_yllg$percent_diff, 3),
      ONT_vs_YLLG_Diffs = ont_vs_yllg$differences,
      ONT_vs_YLLG_Sites = ont_vs_yllg$total_compared,
      ONT_vs_YLLG_Percent = round(ont_vs_yllg$percent_diff, 3),
      Illumina_Coverage = round(illumina_coverage * 100, 1),
      ONT_Coverage = round(ont_coverage * 100, 1),
      stringsAsFactors = FALSE
    ))
  }
}

print(yllg_ref_results, row.names = FALSE)

# TABLE 3: Comparison against Individual 60 Reference (if available)
if(!is.null(individual_60_ref)) {
  cat("\n\n=== TABLE 3: COMPARISON AGAINST INDIVIDUAL 60 REFERENCE ===\n\n")
  
  ind60_ref_results <- data.frame(
    Individual = character(0),
    Illumina_vs_60_Diffs = numeric(0),
    Illumina_vs_60_Sites = numeric(0),
    Illumina_vs_60_Percent = numeric(0),
    ONT_vs_60_Diffs = numeric(0),
    ONT_vs_60_Sites = numeric(0),
    ONT_vs_60_Percent = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for(individual in common_individuals) {
    if(individual %in% names(illumina_seqs) && individual %in% names(ont_seqs)) {
      illumina_seq <- illumina_seqs[[individual]]
      ont_seq <- ont_seqs[[individual]]
      
      # Coverage check
      illumina_chars <- strsplit(illumina_seq, "")[[1]]
      ont_chars <- strsplit(ont_seq, "")[[1]]
      illumina_coverage <- sum(illumina_chars != "N" & illumina_chars != "n") / length(illumina_chars)
      ont_coverage <- sum(ont_chars != "N" & ont_chars != "n") / length(ont_chars)
      
      if(illumina_coverage < 0.6 | ont_coverage < 0.6) next
      
      # Compare to individual 60 reference
      illumina_vs_60 <- compare_to_reference(illumina_seq, individual_60_ref)
      ont_vs_60 <- compare_to_reference(ont_seq, individual_60_ref)
      
      ind60_ref_results <- rbind(ind60_ref_results, data.frame(
        Individual = individual,
        Illumina_vs_60_Diffs = illumina_vs_60$differences,
        Illumina_vs_60_Sites = illumina_vs_60$total_compared,
        Illumina_vs_60_Percent = round(illumina_vs_60$percent_diff, 3),
        ONT_vs_60_Diffs = ont_vs_60$differences,
        ONT_vs_60_Sites = ont_vs_60$total_compared,
        ONT_vs_60_Percent = round(ont_vs_60$percent_diff, 3),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  print(ind60_ref_results, row.names = FALSE)
}

# SUMMARY STATISTICS
cat("\n\n=== SUMMARY STATISTICS ===\n")

cat("\nONT vs Illumina Concordance:\n")
if(nrow(ont_vs_illumina_results) > 0) {
  cat("Mean differences:", round(mean(ont_vs_illumina_results$SNP_Differences, na.rm = TRUE), 2), "SNPs\n")
  cat("Mean percent different:", round(mean(ont_vs_illumina_results$Percent_Different, na.rm = TRUE), 3), "%\n")
  cat("Mean concordance:", round(100 - mean(ont_vs_illumina_results$Percent_Different, na.rm = TRUE), 3), "%\n")
}

cat("\nAccuracy vs Yellow-legged Gull Reference:\n")
if(nrow(yllg_ref_results) > 0) {
  cat("Illumina mean differences:", round(mean(yllg_ref_results$Illumina_vs_YLLG_Diffs, na.rm = TRUE), 2), "SNPs\n")
  cat("ONT mean differences:", round(mean(yllg_ref_results$ONT_vs_YLLG_Diffs, na.rm = TRUE), 2), "SNPs\n")
  cat("Illumina mean accuracy:", round(100 - mean(yllg_ref_results$Illumina_vs_YLLG_Percent, na.rm = TRUE), 3), "%\n")
  cat("ONT mean accuracy:", round(100 - mean(yllg_ref_results$ONT_vs_YLLG_Percent, na.rm = TRUE), 3), "%\n")
  
  # Technology comparison
  illumina_better <- sum(yllg_ref_results$Illumina_vs_YLLG_Diffs < yllg_ref_results$ONT_vs_YLLG_Diffs, na.rm = TRUE)
  ont_better <- sum(yllg_ref_results$ONT_vs_YLLG_Diffs < yllg_ref_results$Illumina_vs_YLLG_Diffs, na.rm = TRUE)
  cat("Individuals where Illumina more accurate:", illumina_better, "\n")
  cat("Individuals where ONT more accurate:", ont_better, "\n")
}


