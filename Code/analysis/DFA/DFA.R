# Main DFA Analysis Script for Management Units
# Simplified - reads existing management unit data instead of recalculating

# Load required libraries
library(dplyr)
library(tidyr)
library(here)
library(sf)
library(ggplot2)
library(RColorBrewer)

############# This bit includes the influence of CPUE, way more fish are caught in the middle, so way more fish from there. 

# Read in the management tributary summary 
mngmt_rivs <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_Units/management_unit_summary.csv")

# Create a wide format timeseries (each row = management river, columns = year-quartile combinations)
timeseries_wide <- mngmt_rivs %>%
  # Create a time period identifier
  mutate(time_period = paste0(year, "_", quartile)) %>%
  # Select key columns
  select(mgmt_river, time_period, percent_of_total_run) %>%
  # Pivot to wide format
  pivot_wider(
    names_from = time_period, 
    values_from = percent_of_total_run,
    values_fill = 0  # Fill missing values with 0
  )

# Display the wide format timeseries
cat("\n=== Wide Format Timeseries ===\n")
print(timeseries_wide)

# Create a long format for easier plotting
timeseries_long <- mngmt_rivs %>%
  mutate(
    # Create decimal year for plotting (Q1=0.125, Q2=0.375, Q3=0.625, Q4=0.875)
    quartile_decimal = case_when(
      quartile == "Q1" ~ 0.125,
      quartile == "Q2" ~ 0.375,
      quartile == "Q3" ~ 0.625,
      quartile == "Q4" ~ 0.875,
      TRUE ~ 0.5
    ),
    year_decimal = year + quartile_decimal
  ) %>%
  arrange(mgmt_river, year_decimal)

# Create a simple line plot showing all management units over time
p1 <- ggplot(timeseries_long, aes(x = year_decimal, y = percent_of_total_run, color = mgmt_river)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(min(timeseries_long$year), max(timeseries_long$year), by = 1),
    labels = scales::number_format(accuracy = 1)
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  labs(
    title = "Management Unit Production Timeseries by Quartile",
    subtitle = "Percent of total watershed run by management river and timing quartile",
    x = "Year (with quarterly increments)",
    y = "Percent of Total Run",
    color = "Management River"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p1)

# Create a faceted plot by management river for clearer individual trends
p2 <- ggplot(timeseries_long, aes(x = year_decimal, y = percent_of_total_run)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 1.5) +
  facet_wrap(~ mgmt_river, scales = "free_y", ncol = 3) +
  scale_x_continuous(
    breaks = seq(min(timeseries_long$year), max(timeseries_long$year), by = 1),
    labels = scales::number_format(accuracy = 1)
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  labs(
    title = "Management Unit Production Timeseries - Individual Rivers",
    subtitle = "Each panel shows one management river across all quarters and years",
    x = "Year",
    y = "Percent of Total Run"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p2)



#################### This portion is the proportion of each quartile made up of each management river. Doesn't take into account that there are different number of fish from each quartile 
# Management Unit Timeseries Analysis
# Shows proportion of each quartile made up of each management river
# Note: Doesn't account for different numbers of fish from each quartile

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Load and prepare data
timeseries_data <- read.csv("Analysis_Results/Management_Units/enhanced_management_unit_analysis.csv")

timeseries_subset <- timeseries_data %>%
  filter(!is.na(mgmt_river), mgmt_river != "", !is.na(proportion_within_quartile)) %>%
  select(year, quartile, mgmt_river, proportion_within_quartile) %>%
  mutate(
    # Create decimal year for smooth plotting (Q1=0.125, Q2=0.375, etc.)
    quartile_decimal = case_when(
      quartile == "Q1" ~ 0.125, quartile == "Q2" ~ 0.375, 
      quartile == "Q3" ~ 0.625, quartile == "Q4" ~ 0.875
    ),
    year_decimal = year + quartile_decimal
  )

# Plot original proportions
num_rivers <- length(unique(timeseries_subset$mgmt_river))

p1 <- ggplot(timeseries_subset, aes(x = year_decimal, y = proportion_within_quartile, color = mgmt_river)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  scale_x_continuous(breaks = seq(min(timeseries_subset$year), max(timeseries_subset$year), by = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = rainbow(num_rivers), name = "Management River") +
  labs(
    title = "Management Unit Proportional Contributions Within Each Quartile",
    subtitle = paste("Relative dominance among", num_rivers, "managed rivers within each time period"),
    x = "Year", y = "Proportion of Managed Rivers"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4))

print(p1)

# Convert to wide format (rivers as rows, time periods as columns)
quartile_timeseries_wide <- timeseries_subset %>%
  mutate(time_period = paste0(year, "_", quartile)) %>%
  select(mgmt_river, time_period, proportion_within_quartile) %>%
  pivot_wider(names_from = time_period, values_from = proportion_within_quartile, values_fill = 0) %>%
  arrange(mgmt_river)

# Z-normalize each river's timeseries (mean=0, sd=1)
river_names <- quartile_timeseries_wide$mgmt_river
numeric_data <- as.matrix(quartile_timeseries_wide[, -1])  # Remove river names column
z_normalized_matrix <- t(scale(t(numeric_data)))           # Z-normalize each row
z_normalized_matrix[is.nan(z_normalized_matrix)] <- 0      # Handle zero-variance rivers

# Combine normalized data with river names and fix column names
quartile_timeseries_z_normalized <- data.frame(
  mgmt_river = river_names,
  z_normalized_matrix,
  stringsAsFactors = FALSE
)

# Fix column names (remove X prefix that R adds to numeric column names)
colnames(quartile_timeseries_z_normalized) <- c("mgmt_river", colnames(quartile_timeseries_wide)[-1])

# Convert z-normalized data back to long format for plotting
z_timeseries_long <- quartile_timeseries_z_normalized %>%
  pivot_longer(cols = -mgmt_river, names_to = "time_period", values_to = "z_score") %>%
  separate(time_period, into = c("year", "quartile"), sep = "_") %>%
  mutate(
    year = as.numeric(year),
    quartile_decimal = case_when(
      quartile == "Q1" ~ 0.125, quartile == "Q2" ~ 0.375, 
      quartile == "Q3" ~ 0.625, quartile == "Q4" ~ 0.875
    ),
    year_decimal = year + quartile_decimal
  )

# Plot z-normalized timeseries
p2 <- ggplot(z_timeseries_long, aes(x = year_decimal, y = z_score, color = mgmt_river)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_x_continuous(breaks = seq(min(z_timeseries_long$year), max(z_timeseries_long$year), by = 1)) +
  scale_y_continuous(
    breaks = seq(-3, 3, by = 1),
    labels = function(x) paste0(x, "σ"),  # Add sigma symbol to show standard deviations
    limits = c(min(z_timeseries_long$z_score) * 1.1, max(z_timeseries_long$z_score) * 1.1)
  ) +
  scale_color_manual(values = rainbow(num_rivers), name = "Management River") +
  labs(
    title = "Z-Normalized Management Unit Timeseries",
    subtitle = "Each river standardized to mean=0, sd=1 • Shows relative timing patterns",
    x = "Year", y = "Z-Score (standard deviations from river's mean)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4))

print(p2)








###################################### 
###################################### Data is ready for DFA Analysis 




