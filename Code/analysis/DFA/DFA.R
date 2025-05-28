# ============================================================================
# Dynamic Factor Analysis (DFA) for Management Unit Salmon Run Timing
# Complete script for both Total Run and Within-Quartile Proportion analysis
# ============================================================================

# Load required libraries --------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(scales)
library(MARSS)
library(gridExtra)
library(grid)  # Added for textGrob function
library(sf)
library(viridis)

# ============================================================================
# PART 1: DATA LOADING AND PREPARATION
# ============================================================================

# Load the long format data
both_ts_long <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")

# Clean up time variables
both_ts_long <- both_ts_long %>%
  mutate(
    Quartile_Clean = case_when(
      grepl("Q1|1", quartile, ignore.case = TRUE) ~ "Q1",
      grepl("Q2|2", quartile, ignore.case = TRUE) ~ "Q2", 
      grepl("Q3|3", quartile, ignore.case = TRUE) ~ "Q3",
      grepl("Q4|4", quartile, ignore.case = TRUE) ~ "Q4",
      TRUE ~ as.character(quartile)
    ),
    Quartile_Num = case_when(
      Quartile_Clean == "Q1" ~ 1,
      Quartile_Clean == "Q2" ~ 2,
      Quartile_Clean == "Q3" ~ 3,
      Quartile_Clean == "Q4" ~ 4,
      TRUE ~ as.numeric(str_extract(quartile, "\\d"))
    ),
    Time_Continuous = year + (Quartile_Num - 1) * 0.25,
    Time_Period = paste0(year, "_", Quartile_Clean)
  )

# Create color palette for rivers
mgmt_rivers <- unique(both_ts_long$mgmt_river)
n_rivers <- length(mgmt_rivers)

if(n_rivers <= 11) {
  river_colors <- RColorBrewer::brewer.pal(max(3, n_rivers), "Spectral")
} else {
  river_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n_rivers)
}
names(river_colors) <- mgmt_rivers

# ============================================================================
# PART 2: VISUALIZATION OF RAW DATA
# ============================================================================

# Plot 1: Total Run Proportion over time
total_run_data <- both_ts_long %>%
  filter(!is.na(total_run_prop), total_run_prop > 0)

p1 <- ggplot(total_run_data, aes(x = Time_Continuous, y = total_run_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(breaks = unique(total_run_data$year)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, NA)) +
  labs(
    title = "Total Run Proportion by Management River Over Time",
    subtitle = "Percentage of entire watershed run from each management river by quartile",
    x = "Year (with quarterly increments)",
    y = "Proportion of Total Watershed Run (%)",
    color = "Management River"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor.x = element_blank()
  ) +
  guides(color = guide_legend(ncol = min(3, ceiling(n_rivers/2))))

# Plot 2: Within-Quartile Proportion over time (FIXED - moved before ggsave)
within_quartile_data <- both_ts_long %>%
  filter(!is.na(within_quartile_prop), within_quartile_prop > 0)

p2 <- ggplot(within_quartile_data, aes(x = Time_Continuous, y = within_quartile_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(breaks = unique(within_quartile_data$year)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, NA)) +
  labs(
    title = "Within-Quartile Proportion by Management River Over Time",
    subtitle = "Proportion of each quartile's production from each management river",
    x = "Year (with quarterly increments)",
    y = "Proportion of Quartile (%)",
    color = "Management River"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor.x = element_blank()
  ) +
  guides(color = guide_legend(ncol = min(3, ceiling(n_rivers/2))))

# Save the initial visualization plots to Figures directory
ggsave(file.path("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures", 
                 "total_run_proportion_timeseries.png"), 
       p1, width = 12, height = 8, dpi = 300, bg = "white")

ggsave(file.path("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures", 
                 "within_quartile_proportion_timeseries.png"), 
       p2, width = 12, height = 8, dpi = 300, bg = "white")

print(p1)
print(p2)

# ============================================================================
# PART 3: PREPARE DATA FOR DFA - TOTAL RUN PROPORTION
# ============================================================================

cat("\n=== PREPARING TOTAL RUN PROPORTION DATA FOR DFA ===\n")

# Convert to wide format for DFA (rivers as rows, time periods as columns)
total_run_wide <- both_ts_long %>%
  filter(!is.na(total_run_prop)) %>%
  select(mgmt_river, Time_Period, total_run_prop) %>%
  pivot_wider(
    names_from = Time_Period, 
    values_from = total_run_prop, 
    values_fill = 0
  )

# Convert to matrix format for MARSS
river_names <- total_run_wide$mgmt_river
total_run_matrix <- as.matrix(total_run_wide[, -1])
rownames(total_run_matrix) <- river_names

# Check for rivers with zero or near-zero production
row_sums <- rowSums(total_run_matrix, na.rm = TRUE)
row_vars <- apply(total_run_matrix, 1, var, na.rm = TRUE)

# Remove rivers with zero production or zero variance
good_rivers <- which(row_sums > 0 & row_vars > 1e-10 & !is.na(row_vars))
total_run_matrix_clean <- total_run_matrix[good_rivers, ]

cat("Total Run Proportion matrix dimensions:", dim(total_run_matrix_clean), "\n")
cat("Rivers included:", rownames(total_run_matrix_clean), "\n")

# ============================================================================
# PART 4: DFA MODEL SELECTION - TOTAL RUN PROPORTION
# ============================================================================

cat("\n=== RUNNING DFA MODEL SELECTION FOR TOTAL RUN PROPORTION ===\n")

aic_results_total <- data.frame(
  states = 1:min(4, nrow(total_run_matrix_clean)),
  AICc = numeric(min(4, nrow(total_run_matrix_clean))),
  stringsAsFactors = FALSE
)

for (m in 1:min(4, nrow(total_run_matrix_clean))) {
  cat("Fitting model with", m, "trends...\n")
  
  model.list <- list(m = m, R = "diagonal and equal")
  
  model <- suppressMessages(suppressWarnings(
    MARSS(total_run_matrix_clean,
          model = model.list,
          z.score = TRUE,
          form = "dfa",
          control = list(maxit = 50000),
          method = "BFGS")
  ))
  
  aic_results_total$AICc[m] <- model$AICc
}

# Calculate delta AICc
aic_results_total$Delta_AICc <- aic_results_total$AICc - min(aic_results_total$AICc, na.rm = TRUE)

cat("\nTotal Run Proportion - Model Selection Results:\n")
print(aic_results_total)

# Find best model
best_m_total <- which.min(aic_results_total$AICc)
cat("\nBest model for Total Run Proportion:", best_m_total, "trends\n")

# ============================================================================
# PART 5: FIT BEST DFA MODEL - TOTAL RUN PROPORTION
# ============================================================================

cat("\n=== FITTING BEST DFA MODEL FOR TOTAL RUN PROPORTION ===\n")

best_model_total <- suppressMessages(suppressWarnings(
  MARSS(total_run_matrix_clean,
        model = list(m = best_m_total, R = "diagonal and equal"),
        z.score = TRUE,
        form = "dfa", 
        control = list(maxit = 50000), 
        method = "BFGS")
))

# Extract trends and loadings
trends_total <- best_model_total$states
Z_total <- coef(best_model_total, type = "matrix")$Z
colnames(Z_total) <- paste("Trend", 1:best_m_total)

# Create time labels for plotting
years <- rep(2017:2021, each = 4)
quarters <- rep(1:4, times = 5)
time_labels <- paste0(years, "-Q", quarters)

# ============================================================================
# PART 6: PLOT DFA RESULTS - TOTAL RUN PROPORTION
# ============================================================================

cat("\n=== PLOTTING DFA RESULTS FOR TOTAL RUN PROPORTION ===\n")

for (i in 1:best_m_total) {
  # Trend time series
  trend_data <- data.frame(
    time = 1:ncol(trends_total),
    value = trends_total[i, ],
    time_label = time_labels[1:ncol(trends_total)]
  )
  
  trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_point(color = "steelblue", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
    labs(title = paste("Trend", i, "Time Series"), x = "Year-Quarter", y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  # Loadings for this trend
  loading_data <- data.frame(Unit = rownames(Z_total), Loading = Z_total[, i]) %>%
    arrange(desc(abs(Loading)))
  
  loadings_plot <- ggplot(loading_data, aes(x = reorder(Unit, Loading), 
                                            y = Loading, fill = Loading > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "coral"),
      name = "", labels = c("Negative", "Positive")
    ) +
    labs(title = paste("Trend", i, "Loadings"), x = "Management Unit", y = "Loading") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12), 
      legend.position = "bottom"
    )
  
  # Display combined plot
  combined_plot <- grid.arrange(trend_plot, loadings_plot, ncol = 2, 
                                top = paste("TOTAL RUN PROPORTION - TREND", i))
  
  # Save the plot to Figures directory
  ggsave(file.path("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures", 
                   paste0("total_run_proportion_trend_", i, "_plot.png")), 
         combined_plot, width = 12, height = 6, dpi = 300, bg = "white")
  
  print(combined_plot)
}

# ============================================================================
# PART 7: PREPARE DATA FOR DFA - WITHIN-QUARTILE PROPORTION
# ============================================================================

cat("\n=== PREPARING WITHIN-QUARTILE PROPORTION DATA FOR DFA ===\n")

# Convert to wide format for DFA (rivers as rows, time periods as columns)
within_quartile_wide <- both_ts_long %>%
  filter(!is.na(within_quartile_prop)) %>%
  select(mgmt_river, Time_Period, within_quartile_prop) %>%
  pivot_wider(
    names_from = Time_Period, 
    values_from = within_quartile_prop, 
    values_fill = 0
  )

# Convert to matrix format for MARSS
river_names_wq <- within_quartile_wide$mgmt_river
within_quartile_matrix <- as.matrix(within_quartile_wide[, -1])
rownames(within_quartile_matrix) <- river_names_wq

# Check for rivers with zero or near-zero production
row_sums_wq <- rowSums(within_quartile_matrix, na.rm = TRUE)
row_vars_wq <- apply(within_quartile_matrix, 1, var, na.rm = TRUE)

# Remove rivers with zero production or zero variance
good_rivers_wq <- which(row_sums_wq > 0 & row_vars_wq > 1e-10 & !is.na(row_vars_wq))
within_quartile_matrix_clean <- within_quartile_matrix[good_rivers_wq, ]

cat("Within-Quartile Proportion matrix dimensions:", dim(within_quartile_matrix_clean), "\n")
cat("Rivers included:", rownames(within_quartile_matrix_clean), "\n")

# ============================================================================
# PART 8: DFA MODEL SELECTION - WITHIN-QUARTILE PROPORTION
# ============================================================================

cat("\n=== RUNNING DFA MODEL SELECTION FOR WITHIN-QUARTILE PROPORTION ===\n")

aic_results_within <- data.frame(
  states = 1:min(4, nrow(within_quartile_matrix_clean)),
  AICc = numeric(min(4, nrow(within_quartile_matrix_clean))),
  stringsAsFactors = FALSE
)

for (m in 1:min(4, nrow(within_quartile_matrix_clean))) {
  cat("Fitting model with", m, "trends...\n")
  
  model.list <- list(m = m, R = "diagonal and equal")
  
  model <- suppressMessages(suppressWarnings(
    MARSS(within_quartile_matrix_clean,
          model = model.list,
          z.score = TRUE,
          form = "dfa",
          control = list(maxit = 50000),
          method = "BFGS")
  ))
  
  aic_results_within$AICc[m] <- model$AICc
}

# Calculate delta AICc
aic_results_within$Delta_AICc <- aic_results_within$AICc - min(aic_results_within$AICc, na.rm = TRUE)

cat("\nWithin-Quartile Proportion - Model Selection Results:\n")
print(aic_results_within)

# Find best model
best_m_within <- which.min(aic_results_within$AICc)
cat("\nBest model for Within-Quartile Proportion:", best_m_within, "trends\n")

# ============================================================================
# PART 9: FIT BEST DFA MODEL - WITHIN-QUARTILE PROPORTION
# ============================================================================

cat("\n=== FITTING BEST DFA MODEL FOR WITHIN-QUARTILE PROPORTION ===\n")

best_model_within <- suppressMessages(suppressWarnings(
  MARSS(within_quartile_matrix_clean,
        model = list(m = best_m_within, R = "diagonal and equal"),
        z.score = TRUE,
        form = "dfa", 
        control = list(maxit = 50000), 
        method = "BFGS")
))

# Extract trends and loadings
trends_within <- best_model_within$states
Z_within <- coef(best_model_within, type = "matrix")$Z
colnames(Z_within) <- paste("Trend", 1:best_m_within)

# ============================================================================
# PART 10: PLOT DFA RESULTS - WITHIN-QUARTILE PROPORTION
# ============================================================================

cat("\n=== PLOTTING DFA RESULTS FOR WITHIN-QUARTILE PROPORTION ===\n")

for (i in 1:best_m_within) {
  # Trend time series
  trend_data_wq <- data.frame(
    time = 1:ncol(trends_within),
    value = trends_within[i, ],
    time_label = time_labels[1:ncol(trends_within)]
  )
  
  trend_plot_wq <- ggplot(trend_data_wq, aes(x = time, y = value)) +
    geom_line(color = "darkgreen", linewidth = 1.2) +
    geom_point(color = "darkgreen", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(breaks = trend_data_wq$time, labels = trend_data_wq$time_label) +
    labs(title = paste("Trend", i, "Time Series"), x = "Year-Quarter", y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  # Loadings for this trend
  loading_data_wq <- data.frame(Unit = rownames(Z_within), Loading = Z_within[, i]) %>%
    arrange(desc(abs(Loading)))
  
  loadings_plot_wq <- ggplot(loading_data_wq, aes(x = reorder(Unit, Loading), 
                                                  y = Loading, fill = Loading > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
      values = c("FALSE" = "darkgreen", "TRUE" = "orange"),
      name = "", labels = c("Negative", "Positive")
    ) +
    labs(title = paste("Trend", i, "Loadings"), x = "Management Unit", y = "Loading") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12), 
      legend.position = "bottom"
    )
  
  # Display combined plot
  combined_plot_wq <- grid.arrange(trend_plot_wq, loadings_plot_wq, ncol = 2, 
                                   top = paste("WITHIN-QUARTILE PROPORTION - TREND", i))
  
  # Save the plot to Figures directory
  ggsave(file.path("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures", 
                   paste0("within_quartile_proportion_trend_", i, "_plot.png")), 
         combined_plot_wq, width = 12, height = 6, dpi = 300, bg = "white")
  
  print(combined_plot_wq)
}

# ============================================================================
# PART 11: SUMMARY AND COMPARISON
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")

# Create summary comparison
summary_df <- data.frame(
  Analysis_Type = c("Total Run Proportion", "Within-Quartile Proportion"),
  Best_Model_Trends = c(best_m_total, best_m_within),
  Best_AICc = c(min(aic_results_total$AICc), min(aic_results_within$AICc)),
  Rivers_Included = c(nrow(total_run_matrix_clean), nrow(within_quartile_matrix_clean)),
  Time_Points = c(ncol(total_run_matrix_clean), ncol(within_quartile_matrix_clean))
)

print(summary_df)

cat("\n=== MODEL SELECTION COMPARISON ===\n")
cat("Total Run Proportion AIC Results:\n")
print(aic_results_total)
cat("\nWithin-Quartile Proportion AIC Results:\n")
print(aic_results_within)

# Save results to files
cat("\n=== SAVING RESULTS ===\n")

# Create output directories
figures_dir <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures"
maps_dir <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Maps"
results_dir <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Results"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(maps_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Save model selection results
write.csv(aic_results_total, file.path(results_dir, "total_run_proportion_model_selection.csv"), row.names = FALSE)
write.csv(aic_results_within, file.path(results_dir, "within_quartile_proportion_model_selection.csv"), row.names = FALSE)

# Save trend data
trends_total_df <- data.frame(
  Time = time_labels[1:ncol(trends_total)],
  as.data.frame(t(trends_total))
)
colnames(trends_total_df)[-1] <- paste("Trend", 1:best_m_total)

trends_within_df <- data.frame(
  Time = time_labels[1:ncol(trends_within)],
  as.data.frame(t(trends_within))
)
colnames(trends_within_df)[-1] <- paste("Trend", 1:best_m_within)

write.csv(trends_total_df, file.path(results_dir, "total_run_proportion_trends.csv"), row.names = FALSE)
write.csv(trends_within_df, file.path(results_dir, "within_quartile_proportion_trends.csv"), row.names = FALSE)

# Save loadings
loadings_total_df <- data.frame(
  Management_Unit = rownames(Z_total),
  Z_total
)
loadings_within_df <- data.frame(
  Management_Unit = rownames(Z_within),
  Z_within
)

write.csv(loadings_total_df, file.path(results_dir, "total_run_proportion_loadings.csv"), row.names = FALSE)
write.csv(loadings_within_df, file.path(results_dir, "within_quartile_proportion_loadings.csv"), row.names = FALSE)

cat("Analysis complete! Results saved to CSV files.\n")

# ============================================================================
# PART 12: CREATE SPATIAL MAPS WITH LOADINGS (FIXED - Stream Order Line Width)
# ============================================================================

cat("\n=== CREATING SPATIAL MAPS WITH LOADINGS ===\n")

# Load spatial data (adjust path as needed)
tryCatch({
  # Load the Kusko edges shapefile with management river information
  edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp", quiet = TRUE)
  basin <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp", quiet = TRUE)
  
  cat("Loaded spatial data successfully\n")
  cat("Edges loaded:", nrow(edges), "features\n")
  cat("Basin loaded:", nrow(basin), "features\n")
  
  # Check if mgmt_river column exists
  if ("mgmt_river" %in% colnames(edges)) {
    cat("Management river data available\n")
    
    # Function to create trend and loading maps
    create_trend_loading_maps <- function(trends_data, Z_matrix, analysis_name, edges_sf, basin_sf) {
      
      n_trends <- ncol(Z_matrix)
      
      for (trend_num in 1:n_trends) {
        cat(paste("Creating maps for", analysis_name, "Trend", trend_num, "\n"))
        
        # Prepare trend time series data
        trend_data <- data.frame(
          time = 1:length(trends_data[trend_num, ]),
          value = as.numeric(trends_data[trend_num, ]),
          time_label = time_labels[1:length(trends_data[trend_num, ])]
        )
        
        # Create trend time series plot
        trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
          geom_line(color = "steelblue", linewidth = 1.5, alpha = 0.8) +
          geom_point(color = "steelblue", size = 3, alpha = 0.9) +
          geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
          scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
          labs(
            title = paste(analysis_name, "- Trend", trend_num, "Time Series"),
            x = "Year-Quarter", 
            y = "Trend Value"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(face = "bold", size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold")
          )
        
        # Prepare loadings data
        loadings_data <- data.frame(
          mgmt_river = rownames(Z_matrix),
          loading = Z_matrix[, trend_num],
          stringsAsFactors = FALSE
        )
        
        # Join loadings with spatial data
        edges_with_loadings <- edges_sf %>%
          left_join(loadings_data, by = "mgmt_river") %>%
          # Only keep edges with management river assignments
          filter(!is.na(mgmt_river) & mgmt_river != "") %>%
          # Calculate line width based on stream order (0.3 to 3.0 range)
          mutate(
            # Assuming stream order column is named 'Str_Order' - adjust if different
            stream_order = ifelse(is.na(Str_Order), 3, Str_Order),  # Default to 3 if missing
            # Scale stream order to line width range (0.3 to 3.0)
            line_width = pmax(0.3, pmin(3.0, 0.3 + (stream_order - min(stream_order, na.rm = TRUE)) * 
                                          (3.0 - 0.3) / (max(stream_order, na.rm = TRUE) - min(stream_order, na.rm = TRUE))))
          )
        
        # Create the spatial map
        if (nrow(edges_with_loadings) > 0) {
          
          # Ensure basin and edges have same CRS
          if (st_crs(basin_sf) != st_crs(edges_with_loadings)) {
            basin_sf <- st_transform(basin_sf, st_crs(edges_with_loadings))
          }
          
          # Get the range of stream orders for legend
          stream_order_range <- range(edges_with_loadings$stream_order, na.rm = TRUE)
          
          # Create the map
          map_plot <- ggplot() +
            # Add basin boundary
            geom_sf(data = basin_sf, fill = "gray95", color = "gray70", 
                    linewidth = 0.5, alpha = 0.3) +
            # Add edges colored by loadings, sized by stream order
            geom_sf(data = edges_with_loadings, 
                    aes(color = loading, linewidth = stream_order), 
                    alpha = 0.8) +
            # Use diverging color scale (blue to red) for loadings
            scale_color_gradient2(
              low = "blue", 
              mid = "white", 
              high = "red", 
              midpoint = 0,
              name = "Loading\nValue",
              limits = c(-max(abs(loadings_data$loading)), max(abs(loadings_data$loading))),
              labels = function(x) sprintf("%.2f", x)
            ) +
            # Scale line width by stream order (0.3 to 3.0 range)
            scale_linewidth_continuous(
              range = c(0.3, 3.0), 
              name = "Stream\nOrder",
              breaks = pretty(stream_order_range, n = 4),
              labels = function(x) as.character(round(x))
            ) +
            coord_sf(datum = NA) +
            labs(
              title = paste(analysis_name, "- Trend", trend_num, "Loadings Map"),
              subtitle = "Rivers colored by loading values (blue = negative, red = positive)",
              caption = "Line thickness represents stream order (higher order = thicker lines)"
            ) +
            theme_void() +
            theme(
              plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
              plot.subtitle = element_text(size = 12, hjust = 0.5),
              plot.caption = element_text(size = 10, hjust = 0.5),
              legend.position = "right",
              legend.title = element_text(face = "bold", size = 10),
              legend.text = element_text(size = 9),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA)
            )
          
          # Save the combined figure to Maps directory (without displaying)
          filename <- file.path("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Maps",
                                paste0(gsub(" ", "_", tolower(analysis_name)), "_trend_", trend_num, "_map.png"))
          
          # Create and save the plot without displaying it
          png(filename, width = 16, height = 8, units = "in", res = 300, bg = "white")
          grid.arrange(
            trend_plot, map_plot, 
            ncol = 2, 
            widths = c(1, 1.2),
            top = textGrob(paste(analysis_name, "- Trend", trend_num, "Analysis"), 
                           gp = gpar(fontsize = 16, fontface = "bold"))
          )
          dev.off()
          cat(paste("Saved:", filename, "\n"))
          
        } else {
          cat(paste("No spatial data available for", analysis_name, "Trend", trend_num, "\n"))
        }
      }
    }
    
    # Create maps for Total Run Proportion analysis
    cat("\nCreating Total Run Proportion maps...\n")
    create_trend_loading_maps(trends_total, Z_total, "Total Run Proportion", edges, basin)
    
    # Create maps for Within-Quartile Proportion analysis  
    cat("\nCreating Within-Quartile Proportion maps...\n")
    create_trend_loading_maps(trends_within, Z_within, "Within-Quartile Proportion", edges, basin)
    
    cat("\nAll spatial maps created successfully!\n")
    
  } else {
    cat("Warning: mgmt_river column not found in spatial data\n")
    cat("Available columns:", paste(colnames(edges), collapse = ", "), "\n")
  }
  
}, error = function(e) {
  cat("Error loading spatial data:", e$message, "\n")
  cat("Spatial maps will be skipped\n")
})

cat("\n=== COMPLETE ANALYSIS WITH SPATIAL MAPS FINISHED ===\n")