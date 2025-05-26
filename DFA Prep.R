# Main DFA Analysis Script
# Simplified and organized DFA calculation

# Load required libraries
library(dplyr)
library(tidyr)
library(here)
library(sf)
library(ggplot2)
library(RColorBrewer)

# Source required files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))
source(here("code/utils/Analysis_utils.R"))

# Parameters
watershed <- "Kusko"
years <- c(2017, 2018, 2019, 2020, 2021)
sensitivity_threshold <- 0.7
min_error <- 0.0006
min_stream_order <- 3
output_dir <- here("Analysis_Results/DFA")

# Main analysis workflow
cat("Starting DFA analysis...\n")

# Process data
results <- process_dfa_data(watershed, years, sensitivity_threshold, min_error, min_stream_order)

# Save files and get z-normalized data
time_series_z_normalized <- save_dfa_files(results, output_dir)

# Create and display plots
plots <- create_dfa_plots(results$time_series_with_total, time_series_z_normalized)

print(plots$overall_timing)
print(plots$huc_contributions)
print(plots$z_normalized_facets)

# Fit DFA model
cat("Fitting DFA model...\n")
dfa_model <- fit_dfa_model(time_series_z_normalized)

cat("Analysis complete!\n")