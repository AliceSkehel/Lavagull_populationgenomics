#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Function to install packages if they don't exist
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    cat("Installing", package, "...\n")
    install.packages(package, dependencies = TRUE, quiet = TRUE)
    library(package, character.only = TRUE)
  } else {
    cat(package, "already installed\n")
  }
}

# Install required packages
cat("=== Setting up R packages for haplotype analysis ===\n")
install_if_missing("ape")
install_if_missing("pegas")
install_if_missing("phangorn")

cat("All packages installed successfully!\n")

