
# mtDNA Platform Comparison Analysis - Minimal Version
library(seqinr)
library(ggplot2)

# Set working directory
setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Main analysis function
analyze_platforms <- function(illumina_file = "illumina_all_angsd.fasta", 
                              ont_file = "ont_all_angsd.fasta") {
  
  # Create standardized combined file
  system(paste0("sed 's/^>LVGU_\\([0-9]*\\)$/>&_Illumina/g' ", illumina_file, " > temp_illumina.fasta"))
  system(paste0("sed 's/^>LVGU_\\([0-9]*\\)$/>&_ONT/g' ", ont_file, " | sed 's/>LVGU_0*\\([0-9]*\\)_ONT/>LVGU_\\1_ONT/g' > temp_ont.fasta"))
  system("cat temp_illumina.fasta temp_ont.fasta > combined.fasta")
  
  # Read sequences
  seqs <- read.fasta("combined.fasta", as.string = TRUE, forceDNAtolower = FALSE)
  
  # Find matched samples
  sample_ids <- gsub("^>|_Illumina$|_ONT$", "", names(seqs))
  matched <- unique(sample_ids[duplicated(sample_ids)])
  
  cat(paste("Found", length(matched), "samples on both platforms\n"))
  
  # Calculate concordance
  results <- data.frame(Sample = matched, Identity = numeric(length(matched)), 
                        N_Illumina = numeric(length(matched)), N_ONT = numeric(length(matched)))
  
  for(i in 1:length(matched)) {
    ill_seq <- seqs[[paste0(">", matched[i], "_Illumina")]]
    ont_seq <- seqs[[paste0(">", matched[i], "_ONT")]]
    
    if(!is.null(ill_seq) && !is.null(ont_seq)) {
      ill_chars <- strsplit(ill_seq, "")[[1]]
      ont_chars <- strsplit(ont_seq, "")[[1]]
      
      results$N_Illumina[i] <- sum(ill_chars == "N")
      results$N_ONT[i] <- sum(ont_chars == "N")
      
      valid_pos <- (ill_chars != "N") & (ont_chars != "N")
      if(sum(valid_pos) > 0) {
        results$Identity[i] <- sum(ill_chars[valid_pos] == ont_chars[valid_pos]) / sum(valid_pos) * 100
      } else {
        results$Identity[i] <- NA
      }
    }
  }
  
  # Summary and plot
  cat(paste("Mean concordance:", round(mean(results$Identity, na.rm = TRUE), 2), "%\n"))
  
  p <- ggplot(results, aes(x = Identity)) + 
    geom_histogram(bins = 15, fill = "skyblue", alpha = 0.7) +
    labs(title = "Platform Concordance", x = "Identity (%)", y = "Count") +
    theme_minimal()
  
    write.csv(results, "concordance_results.csv", row.names = FALSE)
  
  # Cleanup
  system("rm temp_illumina.fasta temp_ont.fasta combined.fasta")
  
  return(results)
}

# Run analysis
results <- analyze_platforms()
print(results)

