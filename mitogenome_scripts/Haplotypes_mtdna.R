# MtDNA Haplotype Analysis
# Alice
# 2025-07-17
# using aligned_v4.fasta

library(ape)
library(pegas)

# Read sequences
sequences <- read.dna("/home/askehel/Sequencing_Combined/mtdna_only/consensus/aligned_v4.fasta", 
                      format = "fasta")

cat("Samples:", nrow(sequences), "\n")
cat("Length:", ncol(sequences), "bp\n\n")

# Create haplotypes - treat N's and gaps as missing data
h <- haplotype(sequences, strict = FALSE)

cat("Number of unique haplotypes:", length(h), "\n\n")

# Get haplotype assignments
hap_labels <- attr(h, "index")
names(hap_labels) <- rownames(sequences)

# Show assignments
cat("Haplotype assignments:\n")
for(i in 1:length(unique(hap_labels))) {
  samples_in_hap <- names(hap_labels)[hap_labels == i]
  cat("Haplotype", i, ":", paste(samples_in_hap, collapse=", "), "\n")
}

# Create network
net <- haploNet(h)

# Plot
pdf("haplotype_network.pdf", width = 10, height = 10)
plot(net, 
     size = attr(net, "freq"), 
     scale.ratio = 2,
     labels = TRUE,
     show.mutation = 2)
title(main = "mtDNA Haplotype Network")
dev.off()

cat("\nNetwork saved: haplotype_network.pdf\n")
