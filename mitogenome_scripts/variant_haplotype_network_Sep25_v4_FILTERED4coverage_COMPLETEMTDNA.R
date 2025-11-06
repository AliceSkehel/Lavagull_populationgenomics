library(seqinr)
library(ape)

# High quality samples based on coverage
good_samples <- c("LVGU_1", "LVGU_2", "LVGU_3", "LVGU_4", "LVGU_5", "LVGU_6",
                  "LVGU_10", "LVGU_11", "LVGU_14", "LVGU_18", "LVGU_20", 
                  "LVGU_23", "LVGU_30", "LVGU_37", "LVGU_40", "LVGU_42", 
                  "LVGU_46", "LVGU_47", "LVGU_48", "LVGU_52", "LVGU_54", 
                  "LVGU_56", "LVGU_60")

# Load sequences
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)
all_seqs <- c(modern_seqs, ancient_seqs)

# Filter to good samples
filtered_seqs <- all_seqs[names(all_seqs) %in% good_samples]

cat("Keeping", length(filtered_seqs), "samples\n")

# Save full mtDNA FASTA
output_dir <- "~/Sequencing_Combined/mtdna_only/mtdna_reextraction/analysis_outputs"
write.fasta(sequences=filtered_seqs, 
            names=names(filtered_seqs),
            file.out=file.path(output_dir, "lava_gull_full_mtdna_coverage_filtered.fasta"))

# Convert to NEXUS
full_dna <- read.dna(file.path(output_dir, "lava_gull_full_mtdna_coverage_filtered.fasta"), 
                     format="fasta")

write.nexus.data(full_dna, 
                 file=file.path(output_dir, "lava_gull_full_mtdna_coverage_filtered.nex"),
                 format="dna")

cat("Saved full mtDNA for", length(filtered_seqs), "coverage-filtered samples\n")

