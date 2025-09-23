# MtDNA Haplotype Analysis
# Alice
# 2025-07-17
# using aligned_v4.fasta

library(ape)

# Read your aligned sequences
sequences <- read.dna("/home/askehel/Sequencing_Combined/mtdna_only/consensus/aligned_v4.fasta", format = "fasta")

# Method 1: Check using distance matrix
dist_matrix <- dist.dna(sequences, model = "raw")
dist_matrix_full <- as.matrix(dist_matrix)

# Find identical sequences (distance = 0)
identical_pairs <- which(dist_matrix_full == 0 & upper.tri(dist_matrix_full), arr.ind = TRUE)

cat("Checking for identical sequences:\n")
if(nrow(identical_pairs) > 0) {
  for(i in 1:nrow(identical_pairs)) {
    cat("Identical:", rownames(dist_matrix_full)[identical_pairs[i,1]], 
        "and", colnames(dist_matrix_full)[identical_pairs[i,2]], "\n")
  }
} else {
  cat("No identical sequences found\n")
}

# Method 2: Convert to strings and compare
seq_strings <- apply(as.character(sequences), 1, paste, collapse = "")
unique_seqs <- unique(seq_strings)

cat("\nTotal samples:", length(seq_strings), "\n")
cat("Unique haplotypes:", length(unique_seqs), "\n")

# Method 3: Check if there are any variable sites
variable_sites <- apply(as.character(sequences), 2, function(x) length(unique(x[x != "-" & x != "n"])) > 1)
num_variable <- sum(variable_sites, na.rm = TRUE)

cat("Number of variable sites:", num_variable, "\n")

# Method 4: Look at the first few positions to see variation
cat("\nFirst 20 positions of each sequence:\n")
for(i in 1:min(5, nrow(sequences))) {
  cat(rownames(sequences)[i], ":", paste(as.character(sequences)[i, 1:20], collapse = ""), "\n")
}

