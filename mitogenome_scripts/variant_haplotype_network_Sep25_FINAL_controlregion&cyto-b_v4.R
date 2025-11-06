library(seqinr)
library(ape)

# Load sequences
modern_seqs <- read.fasta("illumina_all_angsd.fasta", as.string = TRUE)
ancient_seqs <- read.fasta("/home/askehel/ancient_all_angsd.fasta", as.string = TRUE)
all_seqs <- c(modern_seqs, ancient_seqs)

seq_chars <- lapply(all_seqs, function(x) unlist(strsplit(x[1], "")))
seq_length <- length(seq_chars[[1]])

# Filter out problematic samples
bad_samples <- c("LVGU_19", "LVGU_24", "LVGU_476")
filtered_seqs <- all_seqs[!names(all_seqs) %in% bad_samples]
filtered_seq_chars <- seq_chars[!names(seq_chars) %in% bad_samples]

# Also filter by overall quality (keep samples with <30% Ns)
n_percentage <- sapply(filtered_seq_chars, function(seq) {
  sum(seq %in% c("n", "N")) / length(seq) * 100
})
good_samples <- names(filtered_seqs)[n_percentage <= 30]
filtered_seqs <- filtered_seqs[good_samples]
filtered_seq_chars <- filtered_seq_chars[good_samples]

cat("Final sample count:", length(filtered_seqs), "\n")

# Extract control region (positions 15500-16800)
control_start <- 15500
control_end <- 16800

extract_region <- function(seq, start, end) {
  chars <- unlist(strsplit(seq[1], ""))
  paste(chars[start:end], collapse="")
}

control_region_seqs <- lapply(filtered_seqs, extract_region, 
                              start=control_start, end=control_end)
names(control_region_seqs) <- names(filtered_seqs)

# Save outputs
output_dir <- "~/Sequencing_Combined/mtdna_only/mtdna_reextraction/analysis_outputs"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# Save FASTA
write.fasta(sequences=control_region_seqs, 
            names=names(control_region_seqs),
            file.out=file.path(output_dir, "lava_gull_control_region_final.fasta"))

# Convert to NEXUS for PopART
seqs_control <- read.dna(file.path(output_dir, "lava_gull_control_region_final.fasta"), 
                         format="fasta")

write.nexus.data(seqs_control, 
                 file=file.path(output_dir, "lava_gull_control_region_final.nex"),
                 format="dna")

cat("Done. Created control region files (", control_end - control_start + 1, "bp) for", 
    length(control_region_seqs), "samples\n")
cat("NEXUS file for PopART:", file.path(output_dir, "lava_gull_control_region_final.nex"), "\n")
