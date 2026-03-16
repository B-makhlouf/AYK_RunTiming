################################################################################
# VALIDATION: Analysis 2 vs Analysis 4 (Aggregated Across All Years)
################################################################################

library(dplyr)
library(readr)

# File paths
CSV_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Cluster Analysis Results/CSV_Exports"
CLOSURE_END_DOY <- 162  # June 11

# Load data
cumulative_data <- read_csv(file.path(CSV_DIR, "analysis2_cluster_cumulative_by_year.csv"), 
                            show_col_types = FALSE)
closure_data <- read_csv(file.path(CSV_DIR, "analysis4_cluster_closure_protection_by_year.csv"), 
                         show_col_types = FALSE)

# Get cumulative % at closure end (or closest DOY before it) for each year
cumulative_at_closure <- cumulative_data %>%
  filter(doy <= CLOSURE_END_DOY) %>%
  group_by(year, cluster) %>%
  filter(doy == max(doy)) %>%
  ungroup() %>%
  select(year, cluster, cumulative_pct = cluster_cumulative_percent)

# Calculate mean across all years for each cluster
cumulative_avg <- cumulative_at_closure %>%
  group_by(cluster) %>%
  summarise(
    avg_cumulative_pct = mean(cumulative_pct, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

closure_avg <- closure_data %>%
  group_by(cluster) %>%
  summarise(
    avg_closure_pct = mean(closure_protection_pct, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

# Compare
comparison <- cumulative_avg %>%
  left_join(closure_avg, by = "cluster", suffix = c("_cum", "_clo")) %>%
  mutate(
    difference = avg_cumulative_pct - avg_closure_pct,
    match = ifelse(abs(difference) < 0.01, "✓", "✗")
  ) %>%
  select(cluster, avg_cumulative_pct, avg_closure_pct, difference, match, n_years = n_years_cum)

# Print
cat("\n")
cat("Analysis 2 (Cumulative) vs Analysis 4 (Closure Protection)\n")
cat("Average across all years\n")
cat("=" , rep("=", 80), "\n", sep = "")
print(comparison, n = Inf)
cat("=" , rep("=", 80), "\n", sep = "")
cat("\nMax difference:", round(max(abs(comparison$difference)), 4), "%\n")
cat(ifelse(max(abs(comparison$difference)) < 0.01, "✓ PERFECT MATCH\n", "⚠ CHECK DIFFERENCES\n"))
cat("\n")