library(here)
library(tidyverse)

# Load both datasets with original filenames
huc_prod_all <- read.csv(here("/Users/benjaminmakhlouf/Research_repos/03_Shifting-Habitat-Mosaics-II/Analysis_Results/all_huc_production_values.csv"))
runsize <- read.csv(here("/Users/benjaminmakhlouf/Research_repos/03_Shifting-Habitat-Mosaics-II/Data/Escapements_runsize.csv"))

# Process HUC production data for DFA preparation
huc_timeseries <- huc_prod_all %>%
  filter(analysis_type == "DOY") %>%
  mutate(year_quartile = paste0(year, "_", quartile)) %>%
  arrange(Name, year, quartile) %>%
  select(Name, watershed, year, quartile, year_quartile, production_proportion, cpue_percentage)

# Create numeric quartile and time variables
huc_timeseries_numeric <- huc_timeseries %>%
  # Extract numeric value from quartile (1-4)
  mutate(
    quartile_num = as.numeric(gsub("Q", "", quartile)),
    # Create numeric time variable
    time = year + (quartile_num - 1) * 0.25
  )

# Create Q0 rows for each Name and year combination
q0_rows <- huc_timeseries_numeric %>%
  # Get unique Name-year combinations
  select(Name, watershed, year) %>%
  distinct() %>%
  # Create the Q0 rows
  mutate(
    quartile = "Q0",
    quartile_num = 0,
    year_quartile = paste0(year, "_Q0"),
    time = year - 0.25,  # Q0 comes before Q1
    production_proportion = 0,
    cpue_percentage = 0
  )

# Combine original data with Q0 rows
huc_timeseries_complete <- bind_rows(huc_timeseries_numeric, q0_rows) %>%
  arrange(Name, time)

# Plot the time series
ggplot(huc_timeseries_complete, aes(x = time, y = production_proportion, group = Name, color = watershed)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Name, scales = "free_y") +
  theme_bw() +
  # Use year as breaks for cleaner x-axis
  scale_x_continuous(breaks = unique(huc_timeseries_complete$year)) +
  labs(title = "Production Proportion Time Series by HUC",
       x = "Year",
       y = "Production Proportion")

# Combine with runsize data if needed for total run timeseries
# [Add code here to integrate runsize data with your HUC timeseries]


# Combine with runsize data if needed for total run timeseries
# [Add code here to integrate runsize data with your HUC timeseries]







# # Calculate fish numbers per quartile from run size and CPUE data
# cpue_summary <- huc_prod %>%
#   # Remove rows with missing cpue_percentage or quartile
#   filter(!is.na(cpue_percentage), !is.na(quartile)) %>%
#   # Get the first cpue_percentage for each year-watershed-quartile combination
#   group_by(year, watershed, quartile) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   # Join with run size data to get total run for each year
#   left_join(runsize, by = c("year" = "Year", "watershed" = "River")) %>%
#   # Calculate number of fish in each quartile (proportion × total run)
#   mutate(fish_in_quartile = cpue_percentage * Total.Run) %>%
#   select(year, watershed, quartile, cpue_percentage, total_run = Total.Run, fish_in_quartile) %>%
#   arrange(year, quartile)

# # Calculate estimated number of fish from each HUC in each quartile
# huc_fish_estimates <- huc_timeseries %>%
#   # Join with quartile fish numbers
#   left_join(cpue_summary %>% select(year, watershed, quartile, fish_in_quartile), 
#             by = c("year", "watershed", "quartile")) %>%
#   # Calculate estimated fish production per HUC per quartile
#   # (HUC production proportion × fish in that quartile)
#   mutate(estimated_fish_from_huc = production_proportion * fish_in_quartile) %>%
#   # Add helpful metadata
#   mutate(
#     year_quartile = paste0(year, "_", quartile),
#     huc_name = Name
#   ) %>%
#   # Organize columns for clarity
#   select(
#     huc_name = Name,
#     year,
#     watershed,
#     quartile,
#     year_quartile,
#     production_proportion,
#     cpue_percentage,
#     fish_in_quartile,
#     estimated_fish_from_huc
#   ) %>%
#   arrange(huc_name, year, quartile)
# 
# # View the results
# head(huc_fish_estimates, 10)
# 
# # Summary statistics
# cat("Summary of analysis:\n")
# cat("- Number of HUCs:", n_distinct(huc_fish_estimates$huc_name), "\n")
# cat("- Years covered:", min(huc_fish_estimates$year), "to", max(huc_fish_estimates$year), "\n")
# cat("- Watersheds:", paste(unique(huc_fish_estimates$watershed), collapse = ", "), "\n")
# 
# # Verify that quartile fish estimates sum to approximately total run
# quartile_verification <- cpue_summary %>%
#   group_by(year, watershed) %>%
#   summarise(
#     total_run = first(total_run),
#     sum_quartile_fish = sum(fish_in_quartile),
#     difference = total_run - sum_quartile_fish,
#     .groups = "drop"
#   )
# 
# cat("\nVerification - Quartile totals vs Total Run:\n")
# print(quartile_verification)

############ SOMETHINGS WRONG, WILL FIX LATER 

plot_data <- huc_fish_estimates %>%
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
  filter(!is.na(estimated_fish_from_huc))




# 1. Line plot showing all HUCs over time (using decimal year for smooth lines)
p1_continuous <- ggplot(plot_data, aes(x = year_decimal, y = estimated_fish_from_huc, color = huc_name)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_continuous(
    breaks = seq(min(plot_data$year), max(plot_data$year), by = 1),
    labels = scales::number_format(accuracy = 1),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K")) +
  labs(
    title = "Fish Production by HUC Across Quartiles Over Time",
    subtitle = "Estimated number of fish from each HUC by run timing quartile",
    x = "Year (with quarterly increments)",
    y = "Estimated Fish Count (thousands)",
    color = "HUC Name"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank()
  ) +
  guides(color = guide_legend(ncol = 3))


###########



p2_faceted <- ggplot(plot_data, aes(x = factor(year), y = estimated_fish_from_huc, color = huc_name, group = huc_name)) +
  geom_line(size = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(~quartile, scales = "free_y", ncol = 2, 
             labeller = labeller(quartile = c("Q1" = "Q1 (Early Run)", 
                                              "Q2" = "Q2 (Early-Mid Run)", 
                                              "Q3" = "Q3 = (Mid-Late Run)", 
                                              "Q4" = "Q4 (Late Run)"))) +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K")) +
  labs(
    title = "Fish Production by HUC: Seasonal Patterns Over Time",
    subtitle = "Each panel shows production for different run timing quartiles",
    x = "Year",
    y = "Estimated Fish Count (thousands)",
    color = "HUC Name"
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


# 3. Heatmap showing all HUCs, years, and quartiles
p3_heatmap <- plot_data %>%
  ggplot(aes(x = time_period, y = fct_reorder(huc_name, estimated_fish_from_huc, .fun = mean, .desc = TRUE))) +
  geom_tile(aes(fill = estimated_fish_from_huc), color = "white", size = 0.1) +
  scale_fill_viridis_c(
    name = "Fish Count\n(thousands)",
    labels = scales::number_format(scale = 1e-3, suffix = "K"),
    option = "plasma"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Fish Production Heatmap: HUC × Time Period",
    subtitle = "Color intensity represents estimated fish count",
    x = "Time Period (Year-Quartile)",
    y = "HUC (ordered by mean production)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )
