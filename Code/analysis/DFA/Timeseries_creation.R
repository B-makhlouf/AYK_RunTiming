library(here)
library(tidyverse)

# Load management unit data (from your management analysis) and runsize data
mgmt_prod_all <- read.csv(here("Analysis_Results/Management_Units/management_unit_analysis.csv"))
runsize <- read.csv(here("/Users/benjaminmakhlouf/Research_repos/03_Shifting-Habitat-Mosaics-II/Data/Escapements_runsize.csv"))

# Process management unit production data for timeseries analysis
mgmt_timeseries <- mgmt_prod_all %>%
  mutate(year_quartile = paste0(year, "_", quartile)) %>%
  arrange(mgmt_river, year, quartile) %>%
  select(mgmt_river, watershed, year, quartile, year_quartile, percent_of_total_run, total_production, edge_count)

# Create numeric quartile and time variables
mgmt_timeseries_numeric <- mgmt_timeseries %>%
  # Extract numeric value from quartile (1-4)
  mutate(
    quartile_num = as.numeric(gsub("Q", "", quartile)),
    # Create numeric time variable
    time = year + (quartile_num - 1) * 0.25,
    # Convert percent_of_total_run to proportion for consistency
    production_proportion = percent_of_total_run / 100
  )

# Create Q0 rows for each management unit and year combination
q0_rows <- mgmt_timeseries_numeric %>%
  # Get unique mgmt_river-year combinations
  select(mgmt_river, watershed, year) %>%
  distinct() %>%
  # Create the Q0 rows
  mutate(
    quartile = "Q0",
    quartile_num = 0,
    year_quartile = paste0(year, "_Q0"),
    time = year - 0.25,  # Q0 comes before Q1
    percent_of_total_run = 0,
    production_proportion = 0,
    total_production = 0,
    edge_count = 0
  )

# Combine original data with Q0 rows
mgmt_timeseries_complete <- bind_rows(mgmt_timeseries_numeric, q0_rows) %>%
  arrange(mgmt_river, time)

# Plot the time series
ggplot(mgmt_timeseries_complete, aes(x = time, y = production_proportion, group = mgmt_river, color = mgmt_river)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ mgmt_river, scales = "free_y") +
  theme_bw() +
  # Use year as breaks for cleaner x-axis
  scale_x_continuous(breaks = unique(mgmt_timeseries_complete$year)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Management Unit Production Time Series",
       subtitle = "Proportion of total watershed run by management unit over time",
       x = "Year",
       y = "Proportion of Total Run",
       color = "Management River") +
  theme(legend.position = "bottom")

# Calculate fish numbers per quartile from run size and management unit data
# First, get CPUE data by quartile for the watershed
mgmt_summary <- mgmt_timeseries %>%
  # Get total management unit contribution by year-quartile
  group_by(year, watershed, quartile) %>%
  summarise(
    total_mgmt_percent = sum(percent_of_total_run),
    .groups = "drop"
  )

# Join with run size data to calculate absolute fish numbers
mgmt_fish_estimates <- mgmt_timeseries %>%
  # Join with run size data to get total run for each year
  left_join(runsize, by = c("year" = "Year", "watershed" = "River")) %>%
  # Calculate estimated number of fish from each management unit
  # (Management unit % of total run × total run size)
  mutate(
    estimated_fish_from_mgmt = (percent_of_total_run / 100) * Total.Run,
    total_run = Total.Run
  ) %>%
  # Add helpful metadata
  mutate(
    year_quartile = paste0(year, "_", quartile),
    mgmt_name = mgmt_river
  ) %>%
  # Organize columns for clarity
  select(
    mgmt_name = mgmt_river,
    year,
    watershed,
    quartile,
    year_quartile,
    percent_of_total_run,
    total_production,
    edge_count,
    total_run,
    estimated_fish_from_mgmt
  ) %>%
  arrange(mgmt_name, year, quartile)

# View the results
head(mgmt_fish_estimates, 10)

# Summary statistics
cat("Summary of management unit analysis:\n")
cat("- Number of Management Units:", n_distinct(mgmt_fish_estimates$mgmt_name), "\n")
cat("- Years covered:", min(mgmt_fish_estimates$year), "to", max(mgmt_fish_estimates$year), "\n")
cat("- Watersheds:", paste(unique(mgmt_fish_estimates$watershed), collapse = ", "), "\n")

# Verify that management unit estimates are reasonable
mgmt_verification <- mgmt_fish_estimates %>%
  group_by(year, watershed) %>%
  summarise(
    total_run = first(total_run),
    sum_mgmt_fish = sum(estimated_fish_from_mgmt, na.rm = TRUE),
    sum_mgmt_percent = sum(percent_of_total_run, na.rm = TRUE),
    mgmt_coverage = sum_mgmt_fish / total_run * 100,
    .groups = "drop"
  )

cat("\nVerification - Management unit coverage by year:\n")
print(mgmt_verification)

# Create plot data for visualization
plot_data <- mgmt_fish_estimates %>%
  # Create a proper time variable for plotting
  mutate(
    # Create decimal year for continuous time axis (Q1=0.125, Q2=0.375, Q3=0.625, Q4=0.875)
    quartile_numeric = case_when(
      quartile == "Q1" ~ 0.125,
      quartile == "Q2" ~ 0.375,
      quartile == "Q3" ~ 0.625,
      quartile == "Q4" ~ 0.875,
      TRUE ~ 0.5
    ),
    # Decimal year for continuous plotting
    year_decimal = year + quartile_numeric,
    # Factor for discrete plotting
    time_period = fct_inorder(paste0(year, "-", quartile))
  ) %>%
  # Remove any rows with missing fish estimates
  filter(!is.na(estimated_fish_from_mgmt))

# 1. Line plot showing all management units over time (using decimal year for smooth lines)
p1_continuous <- ggplot(plot_data, aes(x = year_decimal, y = estimated_fish_from_mgmt, color = mgmt_name)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_continuous(
    breaks = seq(min(plot_data$year), max(plot_data$year), by = 1),
    labels = scales::number_format(accuracy = 1),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K")) +
  labs(
    title = "Fish Production by Management Unit Across Quartiles Over Time",
    subtitle = "Estimated number of fish from each management unit by run timing quartile",
    x = "Year (with quarterly increments)",
    y = "Estimated Fish Count (thousands)",
    color = "Management River"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank()
  ) +
  guides(color = guide_legend(ncol = 3))

print(p1_continuous)

# 2. Faceted plot by quartile
p2_faceted <- ggplot(plot_data, aes(x = factor(year), y = estimated_fish_from_mgmt, color = mgmt_name, group = mgmt_name)) +
  geom_line(size = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(~quartile, scales = "free_y", ncol = 2, 
             labeller = labeller(quartile = c("Q1" = "Q1 (Early Run)", 
                                              "Q2" = "Q2 (Early-Mid Run)", 
                                              "Q3" = "Q3 (Mid-Late Run)", 
                                              "Q4" = "Q4 (Late Run)"))) +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K")) +
  labs(
    title = "Fish Production by Management Unit: Seasonal Patterns Over Time",
    subtitle = "Each panel shows production for different run timing quartiles",
    x = "Year",
    y = "Estimated Fish Count (thousands)",
    color = "Management River"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(color = guide_legend(ncol = 3))

print(p2_faceted)

# 3. Heatmap showing all management units, years, and quartiles
p3_heatmap <- plot_data %>%
  ggplot(aes(x = time_period, y = fct_reorder(mgmt_name, estimated_fish_from_mgmt, .fun = mean, .desc = TRUE))) +
  geom_tile(aes(fill = estimated_fish_from_mgmt), color = "white", size = 0.1) +
  scale_fill_viridis_c(
    name = "Fish Count\n(thousands)",
    labels = scales::number_format(scale = 1e-3, suffix = "K"),
    option = "plasma"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Fish Production Heatmap: Management Unit × Time Period",
    subtitle = "Color intensity represents estimated fish count",
    x = "Time Period (Year-Quartile)",
    y = "Management River (ordered by mean production)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

print(p3_heatmap)

# 4. Summary plot showing management unit contributions as percentages
p4_percent <- ggplot(plot_data, aes(x = time_period, y = percent_of_total_run, fill = mgmt_name)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_viridis_d(name = "Management\nRiver") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Management Unit Contributions Over Time",
    subtitle = "Stacked percentages of total watershed run by management unit",
    x = "Time Period (Year-Quartile)",
    y = "Percent of Total Run"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 3))

print(p4_percent)

# Save the processed data
write.csv(mgmt_fish_estimates, here("Analysis_Results/Management_Units/management_unit_timeseries_with_fish_estimates.csv"), row.names = FALSE)
cat("\nSaved processed timeseries data with fish estimates\n")