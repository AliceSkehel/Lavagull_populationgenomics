#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Create user library directory if it doesn't exist
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
}

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
  cat("Created user library directory:", user_lib, "\n")
}

# Set library path
.libPaths(c(user_lib, .libPaths()))
cat("Using library path:", .libPaths()[1], "\n")

# Function to install packages in user library
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", package, "in user library...\n")
    install.packages(package, lib = user_lib, dependencies = TRUE)
    library(package, character.only = TRUE)
    cat(package, "installed successfully!\n")
  } else {
    cat(package, "already available\n")
  }
}

# Install required packages
cat("=== Setting up R packages for haplotype analysis ===\n")
install_if_missing("ape")
install_if_missing("pegas") 
install_if_missing("phangorn")

cat("\nAll packages installed successfully in user library!\n")
cat("Library location:", user_lib, "\n")

