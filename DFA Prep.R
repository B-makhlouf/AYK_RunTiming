# Correct DFA Calculation
# Step-by-step calculation of HUC contributions

library(dplyr)
library(tidyr)
library(here)
library(sf)
library(gridExtra)
library(ggpubr)

# Source required files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

# Parameters
watershed <- "Kusko"
years <- c(2017, 2018, 2019, 2020, 2021)
sensitivity_threshold <- 0.7
min_error <- 0.0006
min_stream_order <- 3

# Create output directory
output_dir <- here("Analysis_Results/DFA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize results
all_results <- list()

# Process each year
for (year in years) {
  message(paste("Processing year", year))
  
  # Step 1: Load data
  spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  Huc <- spatial_data$Huc
  natal_data <- load_natal_data(year, watershed)
  
  # Step 2: Calculate total year CPUE
  total_year_cpue <- sum(natal_data$dailyCPUEprop, na.rm = TRUE)
  
  # Step 3: Divide into quartiles
  quartile_data <- divide_doy_quartiles(natal_data)
  quartile_subsets <- quartile_data$subsets
  
  # Step 4: Set up assignment parameters
  pid_iso <- edges$iso_pred
  pid_isose <- edges$isose_pred
  error <- calculate_error(pid_isose, min_error)
  priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
  
  # Step 5: Process each quartile
  for (q in 1:4) {
    message(paste("  Processing Q", q))
    
    current_subset <- quartile_subsets[[q]]
    
    if (nrow(current_subset) == 0) {
      message(paste("    Skipping Q", q, "- no data"))
      next
    }
    
    # Step 5a: Calculate quartile CPUE proportion
    quartile_cpue <- sum(current_subset$dailyCPUEprop, na.rm = TRUE)
    quartile_cpue_proportion <- quartile_cpue / total_year_cpue
    
    # Step 5b: Perform assignment for this quartile
    assignment_matrix <- perform_assignment(
      current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Step 5c: Calculate basin assignments and normalize within quartile
    basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
    total_quartile_assignment <- sum(basin_assign_sum, na.rm = TRUE)
    
    # Proportion each location contributes to THIS quartile (sums to 1.0)
    basin_assign_proportion <- basin_assign_sum / total_quartile_assignment
    
    # Step 5d: Process HUC data
    huc_result <- process_huc_data(edges, basin, Huc, basin_assign_sum, 8)
    
    # Calculate HUC proportions within this quartile (should sum to 1.0)
    total_huc_production <- sum(huc_result$total_production, na.rm = TRUE)
    huc_result$huc_proportion_within_quartile <- huc_result$total_production / total_huc_production
    
    # Step 5e: Calculate final contribution to total year
    # = (HUC proportion within quartile) × (quartile's proportion of total year CPUE)
    huc_result$contribution_to_total_year <- huc_result$huc_proportion_within_quartile * quartile_cpue_proportion
    
    # Step 5f: Store results
    result_data <- st_drop_geometry(huc_result) %>%
      select(Name, huc_proportion_within_quartile, contribution_to_total_year) %>%
      mutate(
        year = year,
        quartile = q,
        quartile_cpue_proportion = quartile_cpue_proportion
      )
    
    key <- paste(year, "Q", q, sep = "_")
    all_results[[key]] <- result_data
    
    message(paste("    Q", q, "CPUE proportion:", round(quartile_cpue_proportion, 3)))
    message(paste("    HUC proportions sum to:", round(sum(huc_result$huc_proportion_within_quartile), 3)))
  }
}

# Step 6: Combine all results
combined_data <- bind_rows(all_results)

# Step 7: Verify calculations
message("\nVerification:")
year_totals <- combined_data %>%
  group_by(year) %>%
  summarise(
    year_cpue_total = sum(quartile_cpue_proportion),
    year_contribution_total = sum(contribution_to_total_year),
    .groups = "drop"
  )
print(year_totals)

quartile_checks <- combined_data %>%
  group_by(year, quartile) %>%
  summarise(
    huc_proportions_sum = sum(huc_proportion_within_quartile),
    .groups = "drop"
  )
print(quartile_checks)

# Step 7b: Load run size data and create absolute run size timeseries
run_size_path <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Escapements_runsize.csv"
if (file.exists(run_size_path)) {
  run_size_data <- read.csv(run_size_path)
  message("\nRun size data loaded successfully")
  
  # Create absolute run size timeseries
  run_size_timeseries <- combined_data %>%
    select(year, quartile, quartile_cpue_proportion) %>%
    distinct() %>%
    left_join(run_size_data %>% select(Year, `Total.Run`), by = c("year" = "Year")) %>%
    mutate(
      absolute_run_size = `Total.Run` * quartile_cpue_proportion,
      time_step = paste0(year, "_Q", quartile)
    ) %>%
    select(time_step, absolute_run_size) %>%
    arrange(time_step)
  
  message("Created absolute run size timeseries")
  
} else {
  warning(paste("Run size file not found at:", run_size_path))
  run_size_timeseries <- NULL
}

# Step 8: Create time series data
time_series <- combined_data %>%
  select(Name, year, quartile, contribution_to_total_year) %>%
  mutate(time_step = paste0(year, "_Q", quartile)) %>%
  select(Name, time_step, contribution_to_total_year)

# Step 8b: Add overall run timeseries as "Total" and absolute run size as "RunSize"
overall_run_timeseries <- combined_data %>%
  select(year, quartile, quartile_cpue_proportion) %>%
  distinct() %>%
  mutate(
    Name = "Total",
    time_step = paste0(year, "_Q", quartile),
    contribution_to_total_year = quartile_cpue_proportion
  ) %>%
  select(Name, time_step, contribution_to_total_year)

# Add absolute run size if available
if (!is.null(run_size_timeseries)) {
  absolute_run_timeseries <- run_size_timeseries %>%
    mutate(
      Name = "RunSize",
      contribution_to_total_year = absolute_run_size / 1000  # Convert to thousands for better scale
    ) %>%
    select(Name, time_step, contribution_to_total_year)
  
  # Combine HUC data with Total and RunSize
  time_series_with_total <- bind_rows(time_series, overall_run_timeseries, absolute_run_timeseries)
} else {
  # Just combine HUC data with Total
  time_series_with_total <- bind_rows(time_series, overall_run_timeseries)
}

# Step 9: Save data
long_filepath <- file.path(output_dir, "kusko_correct_dfa_long.csv")
write.csv(combined_data, long_filepath, row.names = FALSE)

wide_filepath <- file.path(output_dir, "kusko_correct_dfa_wide.csv")
time_series_wide <- time_series_with_total %>%
  pivot_wider(names_from = time_step, values_from = contribution_to_total_year, values_fill = 0)
write.csv(time_series_wide, wide_filepath, row.names = FALSE)

message(paste("\nSaved corrected data to:"))
message(paste("  Long format:", long_filepath))
message(paste("  Wide format (includes Total):", wide_filepath))

# Step 10: Create plots
plot_data_all <- time_series_with_total %>%
  separate(time_step, into = c("year", "quartile"), sep = "_Q", remove = FALSE) %>%
  mutate(
    year_numeric = as.numeric(year),
    quartile_numeric = as.numeric(quartile),
    time_continuous = year_numeric + (quartile_numeric - 1) * 0.25
  )

# Plot 1: Overall run timing (Total only)
plot_data_total <- plot_data_all %>% filter(Name == "Total")

p1 <- ggplot(plot_data_total, aes(x = time_continuous, y = contribution_to_total_year * 100)) +
  geom_line(linewidth = 2, color = "black") +
  geom_point(size = 3, color = "black") +
  scale_x_continuous(breaks = seq(2017, 2021, 1), labels = seq(2017, 2021, 1)) +
  labs(
    title = "Overall Salmon Run Timing Pattern",
    subtitle = "Proportion of annual CPUE by quartile",
    x = "Year (with quarterly increments)",
    y = "Proportion of Annual Run (%)",
    caption = "Shows seasonal timing pattern across years"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    plot.margin = margin(5, 5, 5, 5, "mm")  # Add consistent margins
  )

# Plot 2: Absolute run size (RunSize only, if available)
if (!is.null(run_size_timeseries)) {
  plot_data_runsize <- plot_data_all %>% filter(Name == "RunSize")
  
  p2 <- ggplot(plot_data_runsize, aes(x = time_continuous, y = contribution_to_total_year)) +
    geom_line(linewidth = 2, color = "blue") +
    geom_point(size = 3, color = "blue") +
    scale_x_continuous(breaks = seq(2017, 2021, 1), labels = seq(2017, 2021, 1)) +
    labs(
      title = "Absolute Run Size by Quartile",
      subtitle = "Number of fish (in thousands) by quartile",
      x = "Year (with quarterly increments)",
      y = "Run Size (thousands of fish)",
      caption = "Total run × quartile CPUE proportion"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
} else {
  p2 <- NULL
}

# Plot 3: Individual HUC contributions (no legend, aligned with other plots)
plot_data_hucs <- plot_data_all %>% filter(Name != "Total" & Name != "RunSize")

p3 <- ggplot(plot_data_hucs, aes(x = time_continuous, y = contribution_to_total_year * 100, color = Name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(plot_data_hucs$Name))), 
                     name = "HUC") +
  scale_x_continuous(breaks = seq(2017, 2021, 1), labels = seq(2017, 2021, 1)) +
  labs(
    title = "Individual HUC Contributions to Annual Salmon Run",
    subtitle = "Each HUC's contribution by quartile",
    x = "Year (with quarterly increments)",
    y = "Contribution to Annual Run (%)",
    caption = "All HUCs combined should sum to 100% for each year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "none",  # Remove legend from this plot
    plot.margin = margin(5, 5, 5, 5, "mm")  # Same margins as other plots
  )

# Extract legend from p3 before removing it
p3_with_legend <- ggplot(plot_data_hucs, aes(x = time_continuous, y = contribution_to_total_year * 100, color = Name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(plot_data_hucs$Name))), 
                     name = "HUC") +
  theme_minimal() +
  theme(legend.position = "right")

legend <- get_legend(p3_with_legend)

# Combine plots with shared legend
if (!is.null(p2)) {
  combined_plot <- ggarrange(
    ggarrange(p1, p2, p3, ncol = 1, heights = c(1, 1, 1), align = "v"),
    legend,
    ncol = 2,
    widths = c(4, 1)
  )
} else {
  combined_plot <- ggarrange(
    ggarrange(p1, p3, ncol = 1, heights = c(1, 1), align = "v"),
    legend,
    ncol = 2,
    widths = c(4, 1)
  )
}

# Save plots
figures_dir <- here("Figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Save individual plots
plot1_path <- file.path(figures_dir, "overall_run_timing.pdf")
plot3_path <- file.path(figures_dir, "huc_contributions.pdf")
combined_path <- file.path(figures_dir, "combined_run_analysis.pdf")

ggsave(plot1_path, p1, width = 12, height = 6, device = "pdf")
ggsave(plot3_path, p3, width = 14, height = 8, device = "pdf")
ggsave(combined_path, combined_plot, width = 14, height = ifelse(!is.null(p2), 15, 12), device = "pdf")

if (!is.null(p2)) {
  plot2_path <- file.path(figures_dir, "absolute_run_size.pdf")
  ggsave(plot2_path, p2, width = 12, height = 6, device = "pdf")
  
  message(paste("Saved plots to:"))
  message(paste("  Overall run timing:", plot1_path))
  message(paste("  Absolute run size:", plot2_path))
  message(paste("  HUC contributions:", plot3_path))
  message(paste("  Combined (stacked):", combined_path))
} else {
  message(paste("Saved plots to:"))
  message(paste("  Overall run timing:", plot1_path))
  message(paste("  HUC contributions:", plot3_path))
  message(paste("  Combined (stacked):", combined_path))
  message("  Note: Absolute run size plot not created (run size data not found)")
}



##################################################
##################################################

# DFA Analysis Setup
# Clear environment and load data

rm(list = ls())

library(here)
library(dplyr)

# Load the timeseries data for the HUCs
huc_timeseries <- read.csv(here("Analysis_Results", "DFA", "kusko_correct_dfa_wide.csv"), stringsAsFactors = FALSE)

# remove the "Middle Fork Kuskokwim River", which is 0 

huc_timeseries <- huc_timeseries %>%
  filter(!grepl("Middle Fork Kuskokwim River", Name))

# Extract the first column (names) and the rest (ts values)
huc_names <- huc_timeseries[, 1]
ts_values <- huc_timeseries[, -1]

# Z-normalize each ROW (each HUC time series)
z_normalized_ts <- t(scale(t(ts_values)))

# Check that each row is normalized
rowMeans(z_normalized_ts, na.rm = TRUE)  # Should be ~0 for each HUC
apply(z_normalized_ts, 1, sd, na.rm = TRUE)  # Should be ~1 for each HUC


dfa_data <- t(z_normalized_ts)


n_trends <- 3  # Number of common trends to extract
n_ts <- nrow(z_normalized_ts)  # Number of time series (HUCs)

# Step 4: Define the DFA model structure
# Z matrix (loadings) - how each time series loads on the common trends
Z_model <- matrix(list(0), n_ts, n_trends)
Z_model[1:n_ts, 1] <- paste0("z", 1:n_ts)  # Each HUC gets its own loading

# R matrix (observation errors) - assume independent errors
R_model <- "diagonal and unequal"

# Step 5: Create the model list
model_list <- list(
  Z = Z_model,
  A = "scaling",  # No level effects
  R = R_model,
  B = "identity", # Random walk for trends
  U = "zero",     # No trend in the trends
  Q = "identity", # Process variance for trends
  x0 = "zero"     # Initial state
)

library(MARSS)
cat("Fitting DFA with", n_trends, "common trend(s)...\n")
dfa_fit <- MARSS(dfa_data, model = model_list, 
                 control = list(maxit = 500))
