# log likelihoods of ngsadmix with the k= 10 repeated x6 times
# Askehel
# 4 th November 25'

library(tidyverse)

# Set working directory
setwd("~/Sequencing_Combined/angsd_variant_calling_all35/diversity_analysis/")

# Read the summary data
summary_data <- read.table("ngsadmix_summary.txt", header=TRUE)

# Calculate SD for each K from the full likelihood data
likelihood_data <- read.table("ngsadmix_likelihoods.txt", header=TRUE)

# Calculate mean and SD for each K
k_stats <- likelihood_data %>%
  group_by(K) %>%
  summarise(
    mean_like = mean(LogLikelihood),
    sd_like = sd(LogLikelihood),
    best_like = max(LogLikelihood)
  )


plot(k_stats$K, k_stats$mean_like/1e6, 
     type="b", 
     pch=19,
     xlab="K",
     ylab="Log Likelihood (×10⁶)",
     main="NGSadmix Likelihood by K",
     ylim=c(min(k_stats$mean_like - k_stats$sd_like)/1e6, 
            max(k_stats$mean_like + k_stats$sd_like)/1e6),
     las=1,
     xaxt="n")  # Suppress default x-axis

# Add custom x-axis with every integer
axis(1, at=1:10, labels=1:10)

# Add error bars (mean ± SD)
arrows(k_stats$K, 
       (k_stats$mean_like - k_stats$sd_like)/1e6,
       k_stats$K,
       (k_stats$mean_like + k_stats$sd_like)/1e6,
       angle=90, code=3, length=0.1)

grid()

# Check the variation for K=1 and K=2
likelihood_data %>%
  filter(K %in% c(1, 2)) %>%
  group_by(K) %>%
  summarise(
    mean = mean(LogLikelihood),
    sd = sd(LogLikelihood),
    min = min(LogLikelihood),
    max = max(LogLikelihood)
  )


