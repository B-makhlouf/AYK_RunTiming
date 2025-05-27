# Main DFA Analysis Script for Management Units
# Simplified - reads existing management unit data instead of recalculating

# Load required libraries
library(dplyr)
library(tidyr)
library(here)
library(sf)
library(ggplot2)
library(RColorBrewer)


###################### DFA of proportion * CPUE ######################
#######################################################################
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

# make timeseries_wide into a matrix 
river_names<- timeseries_wide$mgmt_river  # Extract river names
timeseries_matrix <- as.matrix(timeseries_wide[,-1])  # Exclude the first column (mgmt_river)
rownames(timeseries_matrix) <- river_names  # Set row names to river names

state_numbers <- 1:10
aicc_values <- numeric(length(state_numbers))

for (i in seq_along(state_numbers)) {
  m <- state_numbers[i]
  model.list <- list(m = m, R = "diagonal and equal")
  
  model <- suppressMessages(suppressWarnings(
    MARSS(timeseries_matrix,  # Using your existing matrix
          model = model.list,
          z.score = FALSE,
          form = "dfa", 
          control = list(maxit = 50000), 
          method = "BFGS")
  ))
  
  aicc_values[i] <- model$AICc
}

# Results
results <- data.frame(
  states = state_numbers,
  AICc = aicc_values,
  delta_AICc = aicc_values - min(aicc_values)
)

print(results)

# Fit the optimal 7-trend model
best_model <- suppressMessages(suppressWarnings(
  MARSS(timeseries_matrix,
        model = list(m = 6, R = "diagonal and equal"),
        z.score = FALSE,
        form = "dfa", 
        control = list(maxit = 50000), 
        method = "BFGS")
))

# Extract trends and loadings
trends <- best_model$states
Z <- coef(best_model, type = "matrix")$Z
colnames(Z) <- paste("Trend", 1:6)

# Create plots for each trend (side by side)
plot_list <- list()

# Create time labels - assuming quarterly data from 2017-2021
years <- rep(2017:2021, each = 4)
quarters <- rep(1:4, times = 5)
time_labels <- paste0(years, "-Q", quarters)

for (i in 1:6) {
  # Trend time series with proper time labels
  trend_data <- data.frame(
    time = 1:ncol(trends),
    year = years[1:ncol(trends)],
    quarter = quarters[1:ncol(trends)],
    value = trends[i, ],
    time_label = time_labels[1:ncol(trends)]
  )
  
  trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_point(color = "steelblue", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(
      breaks = trend_data$time,
      labels = trend_data$time_label
    ) +
    labs(title = paste("Trend", i, "Time Series"),
         x = "Year-Quarter", y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  # Loadings for this trend
  loading_data <- data.frame(
    Unit = rownames(Z),
    Loading = Z[, i]
  ) %>%
    arrange(desc(abs(Loading)))
  
  loadings_plot <- ggplot(loading_data, aes(x = reorder(Unit, Loading), y = Loading, 
                                            fill = Loading > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "coral"),
                      name = "", labels = c("Negative", "Positive")) +
    labs(title = paste("Trend", i, "Loadings"),
         x = "Management Unit", y = "Loading") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom")
  
  # Combine trend and loadings side by side
  combined_plot <- grid.arrange(trend_plot, loadings_plot, ncol = 2, 
                                top = paste("TREND", i))
  
  plot_list[[i]] <- combined_plot
}

# Display all trends (you can also save them individually)
for (i in 1:6) {
  print(plot_list[[i]])
}

# Also create a summary table showing the strongest loadings for each trend
cat("\n=== SUMMARY: STRONGEST LOADINGS FOR EACH TREND ===\n")
for (i in 1:6) {
  loadings_for_trend <- Z[, i]
  
  # Sort by absolute value to find strongest relationships
  sorted_loadings <- sort(loadings_for_trend, decreasing = TRUE)
  
  cat(paste0("\nTrend ", i, ":\n"))
  cat("  Most POSITIVE: ")
  pos_units <- names(sorted_loadings[sorted_loadings > 0])[1:3]
  pos_values <- round(sorted_loadings[sorted_loadings > 0][1:3], 2)
  cat(paste(paste0(pos_units, " (", pos_values, ")"), collapse = ", "), "\n")
  
  cat("  Most NEGATIVE: ")
  neg_loadings <- sort(loadings_for_trend[loadings_for_trend < 0])
  if(length(neg_loadings) > 0) {
    neg_units <- names(neg_loadings)[1:min(3, length(neg_loadings))]
    neg_values <- round(neg_loadings[1:min(3, length(neg_loadings))], 2)
    cat(paste(paste0(neg_units, " (", neg_values, ")"), collapse = ", "), "\n")
  } else {
    cat("None\n")
  }
}



















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

# Convert to a matrix
river_names <- quartile_timeseries_z_normalized$mgmt_river
dfa_matrix <- as.matrix(quartile_timeseries_z_normalized[, -1])  # Remove river names column
rownames(dfa_matrix) <- river_names

###################################### 
###################################### Data is ready for DFA Analysis 

state_numbers <- 1:10
aicc_values <- numeric(length(state_numbers))

for (i in seq_along(state_numbers)) {
  m <- state_numbers[i]
  model.list <- list(m = m, R = "diagonal and equal")
  
  model <- suppressMessages(suppressWarnings(
    MARSS(dfa_matrix,  # Using your existing matrix
          model = model.list,
          z.score = FALSE,
          form = "dfa", 
          control = list(maxit = 50000), 
          method = "BFGS")
  ))
  
  aicc_values[i] <- model$AICc
}

# Results
results <- data.frame(
  states = state_numbers,
  AICc = aicc_values,
  delta_AICc = aicc_values - min(aicc_values)
)

print(results)

# DFA Analysis: Trends and Loadings Side by Side
library(MARSS)
library(tidyverse)
library(ggplot2)
library(gridExtra)

# Fit the optimal 6-trend model
best_model <- suppressMessages(suppressWarnings(
  MARSS(dfa_matrix,
        model = list(m = 4, R = "diagonal and equal"),
        z.score = FALSE,
        form = "dfa", 
        control = list(maxit = 50000), 
        method = "BFGS")
))

# Extract trends and loadings
trends <- best_model$states
Z <- coef(best_model, type = "matrix")$Z
colnames(Z) <- paste("Trend", 1:4)

# Create plots for each trend (side by side)
plot_list <- list()

# Create time labels - assuming quarterly data from 2017-2021
years <- rep(2017:2021, each = 4)
quarters <- rep(1:4, times = 5)
time_labels <- paste0(years, "-Q", quarters)

for (i in 1:4) {
  # Trend time series with proper time labels
  trend_data <- data.frame(
    time = 1:ncol(trends),
    year = years[1:ncol(trends)],
    quarter = quarters[1:ncol(trends)],
    value = trends[i, ],
    time_label = time_labels[1:ncol(trends)]
  )
  
  trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_point(color = "steelblue", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(
      breaks = trend_data$time,
      labels = trend_data$time_label
    ) +
    labs(title = paste("Trend", i, "Time Series"),
         x = "Year-Quarter", y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  # Loadings for this trend
  loading_data <- data.frame(
    Unit = rownames(Z),
    Loading = Z[, i]
  ) %>%
    arrange(desc(abs(Loading)))
  
  loadings_plot <- ggplot(loading_data, aes(x = reorder(Unit, Loading), y = Loading, 
                                            fill = Loading > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "coral"),
                      name = "", labels = c("Negative", "Positive")) +
    labs(title = paste("Trend", i, "Loadings"),
         x = "Management Unit", y = "Loading") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom")
  
  # Combine trend and loadings side by side
  combined_plot <- grid.arrange(trend_plot, loadings_plot, ncol = 2, 
                                top = paste("TREND", i))
  
  plot_list[[i]] <- combined_plot
}

# Display all trends (you can also save them individually)
for (i in 1:6) {
  print(plot_list[[i]])
}
