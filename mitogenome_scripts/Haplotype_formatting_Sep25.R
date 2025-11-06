## Askehel
# 23rd September 2025
# Trying again to make a haplotype network of my mtDNA
# I am going to do most of the tree building using an iqtree in Figtree
# First I need to format the illumina, ancient and ont mtDNA to be ready to go into the command line
# Then I need to pull the file out onto my desktop where I'll then open it into FigTree

# PREP YOUR mtDNA DATA FOR IQ-TREE
library(seqinr)

# Set working directory 
setwd("/home/askehel/Sequencing_Combined/mtdna_only/mtdna_reextraction/")

# Load your sequences (your existing code)
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)

# Combine sequences
seqs <- c(modern_seqs, ancient_seqs)

# Clean up sample names 
names(seqs) <- gsub("LVGU_", "", names(seqs))

# Filter by coverage (your 60% threshold)
coverage <- sapply(seqs, function(x) {
  chars <- strsplit(x, "")[[1]]
  sum(chars != "N" & chars != "n") / length(chars)
})
seqs <- seqs[coverage >= 0.6]

# EXPORT FOR IQ-TREE (save in analysis_outputs directory)
write.fasta(sequences = seqs, 
            names = names(seqs), 
            file.out = "analysis_outputs/mtdna_for_iqtree.fasta")

# Print summary
cat("Data prepared for IQ-TREE!\n")
cat("File created: analysis_outputs/mtdna_for_iqtree.fasta\n")
cat("Number of sequences:", length(seqs), "\n")
cat("Sequence length:", nchar(seqs[[1]]), "bp\n")
cat("Coverage threshold: 60%\n")

cat("\nNext step: Run in terminal:\n")
cat("cd", file.path(getwd(), "analysis_outputs"), "\n")
cat("iqtree -s mtdna_for_iqtree.fasta -m MFP -bb 1000 -nt AUTO\n")

