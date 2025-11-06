
pacman::p_load(ggplot2,dplyr,tidyr,RColorBrewer)

# Load libraries

# Set working directory to where the DATA is (outside git)
setwd("~/Sequencing_Combined/angsd_variant_calling_all35/ngsadmix_analysis")

# Load metadata
metadata <- read.csv("~/Sequencing_Combined/angsd_variant_calling_all35/Kinship_Summary_Data.csv")

# Check how many samples
cat("Metadata has", nrow(metadata), "samples\n")

cluster_colors <- RColorBrewer::brewer.pal(6, "Pastel1")

# Read K=2 to check dimensions
k2 <- read.table("lava_gulls_k2.qopt")
cat("K=2 qopt has", nrow(k2), "samples\n")

# Function to create admixture barplot
plot_k <- function(k_val, colors = NULL) {
  # Read admixture proportions
  q <- read.table(paste0("lava_gulls_k", k_val, ".qopt"))
  
  # Add sample info
  q$Sample <- metadata$ID[1:nrow(q)]
  q$Island <- metadata$ISLAND.or.Country[1:nrow(q)]
  
  # Order by island then by admixture
  q <- q %>% arrange(Island, V1)
  q$Order <- 1:nrow(q)
  
  # Reshape for plotting
  q_long <- q %>%
    select(Sample, Island, Order, starts_with("V")) %>%
    pivot_longer(cols = starts_with("V"), 
                 names_to = "Cluster", 
                 values_to = "Proportion")
  
  # Set default colors if not provided
  if(is.null(colors)) {
    colors <- cluster_colors[1:k_val]
  }
  
  # Plot
  p <- ggplot(q_long, aes(x = Order, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = colors) +
    facet_grid(. ~ Island, scales = "free_x", space = "free_x") +
    theme_minimal() +
    labs(title = paste0("K = ", k_val),
         x = "Individuals",
         y = "Ancestry Proportion") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          strip.text = element_text(size = 10, angle = 90),
          legend.position = "none")
  
  return(p)
}
# Create plots
k2_plot <- plot_k(2)
k3_plot <- plot_k(3)
k4_plot <- plot_k(4)
k5_plot <- plot_k(5)
k6_plot <- plot_k(6)

# Print to screen
print(k2_plot)
print(k3_plot)
print(k4_plot)
print(k5_plot)
print(k6_plot)

# Example: a CSV with columns K, Run, LogLikelihood
log_dir <- "~/Sequencing_Combined/angsd_variant_calling_all35/ngsadmix_analysis"

# List all log files
log_files <- list.files(log_dir, pattern = "lava_gulls_k[0-9]+\\.log$", full.names = TRUE)

# Function to extract K and log-likelihood
extract_likelihood <- function(file) {
  # Read log file lines
  lines <- readLines(file)
  
  # Extract K value from filename (e.g. k3)
  k_val <- as.numeric(str_extract(basename(file), "(?<=k)\\d+"))
  
  # Find line containing 'best like=' or 'like='
  like_line <- lines[str_detect(lines, "best like|like=")]
  
  # Extract numeric likelihood value
  loglik <- as.numeric(str_extract(like_line, "-?\\d+\\.?\\d*"))
  
  # Return data frame
  data.frame(K = k_val, LogLikelihood = loglik)
}

# Apply to all files
likelihood_df <- map_dfr(log_files, extract_likelihood)

# Check results
print(likelihood_df)


# Summarise mean and SD per K
likelihood_summary <- likelihood_df %>%
  group_by(K) %>%
  summarise(
    mean_loglik = mean(LogLikelihood),
    sd_loglik = sd(LogLikelihood),
    .groups = "drop"
  )

# Plot mean ± SD for each K
likelihood_plot <- ggplot(likelihood_summary, aes(x = K, y = mean_loglik)) +
  geom_line(linewidth = 1, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = mean_loglik - sd_loglik,
                    ymax = mean_loglik + sd_loglik),
                width = 0.15, color = "#0072B2") +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of ancestral populations (K)",
    y = "Mean log likelihood (± SD)",
    title = "Likelihood values from NGSadmix"
  ) +
  scale_x_continuous(breaks = unique(likelihood_summary$K)) +
  scale_y_reverse()  # Optional, flips axis so lower = better fit

# Print the plot
print(likelihood_plot)


# Your likelihood data
likelihood_df <- data.frame(
  K = 2:6,
  LogLikelihood = c(-20276549, -19597709, -19090301, -18663480, -18252468)
)

# Compute ΔLogLikelihood between consecutive K
likelihood_df <- likelihood_df %>%
  mutate(Delta = c(NA, diff(LogLikelihood)))

# Scree plot
ggplot(likelihood_df, aes(x = K, y = LogLikelihood)) +
  geom_line(color = "#0072B2", linewidth = 1) +
  geom_point(size = 3, color = "#0072B2") +
  geom_text(aes(label = round(Delta, 0)), vjust = -1.2, color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of clusters (K)",
    y = "Log-likelihood",
    title = "NGSadmix likelihood Scree plot"
  ) +
  scale_x_continuous(breaks = 2:6) +
  scale_y_reverse()  # flips axis so lower = better fit

     