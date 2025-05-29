# ============================================================================
# Dynamic Factor Analysis (DFA) for Salmon Proportion of Total
# Polished Analysis Script - Model comparison and visualization
# ============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(MARSS)
  library(gridExtra)
  library(sf)
  library(grid)
})

# Set paths
DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
FIGURE_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures"
MAP_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Maps"
SPATIAL_DATA_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"

# ============================================================================
# 1. DATA PREPARATION
# ============================================================================

cat("Loading and preparing data...\n")

# Load and clean data
both_ts_long <- read.csv(DATA_PATH) %>%
  mutate(
    Quartile_Clean = case_when(
      grepl("Q1|1", quartile, ignore.case = TRUE) ~ "Q1",
      grepl("Q2|2", quartile, ignore.case = TRUE) ~ "Q2", 
      grepl("Q3|3", quartile, ignore.case = TRUE) ~ "Q3",
      grepl("Q4|4", quartile, ignore.case = TRUE) ~ "Q4",
      TRUE ~ as.character(quartile)
    ),
    Quartile_Num = case_when(
      Quartile_Clean == "Q1" ~ 1, Quartile_Clean == "Q2" ~ 2,
      Quartile_Clean == "Q3" ~ 3, Quartile_Clean == "Q4" ~ 4,
      TRUE ~ as.numeric(str_extract(quartile, "\\d"))
    ),
    Time_Continuous = year + (Quartile_Num - 1) * 0.25,
    Time_Period = paste0(year, "_", Quartile_Clean)
  )

# Add Q0 entries (set to 0 before all Q1s)
q0_entries <- both_ts_long %>%
  select(mgmt_river, year) %>%
  distinct() %>%
  mutate(
    quartile = "Q0",
    Quartile_Clean = "Q0",
    Quartile_Num = 0,
    Time_Continuous = year - 0.25,
    Time_Period = paste0(year, "_Q0"),
    total_run_prop = 0
  )

# Add missing columns with appropriate defaults
missing_cols <- setdiff(names(both_ts_long), names(q0_entries))
for(col in missing_cols) {
  if(is.numeric(both_ts_long[[col]])) {
    q0_entries[[col]] <- 0
  } else {
    q0_entries[[col]] <- NA
  }
}

# Combine data
both_ts_long_with_q0 <- bind_rows(both_ts_long, q0_entries) %>%
  arrange(mgmt_river, year, Quartile_Num)

cat(sprintf("Original data: %d rows | With Q0: %d rows | Added: %d Q0 entries\n", 
            nrow(both_ts_long), nrow(both_ts_long_with_q0), nrow(q0_entries)))

# Filter valid data
proportion_total_data <- both_ts_long_with_q0 %>%
  filter(!is.na(total_run_prop))

# Define river colors
mgmt_rivers <- unique(proportion_total_data$mgmt_river)
n_rivers <- length(mgmt_rivers)
river_colors <- if (n_rivers <= 11) {
  brewer.pal(max(3, n_rivers), "Spectral")
} else {
  colorRampPalette(brewer.pal(11, "Spectral"))(n_rivers)
}
names(river_colors) <- mgmt_rivers

# Plot
p1 <- ggplot(proportion_total_data, aes(x = Time_Continuous, y = total_run_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(breaks = unique(proportion_total_data$year)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, NA)) +
  labs(
    title = "Proportion of Total by Management River Over Time",
    x = "Year (Quarterly)",
    y = "Proportion of Total (%)",
    color = "Management River"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# Save and print
ggsave(file.path(FIGURE_PATH, "total_run_prop_timeseries.png"), 
       p1, width = 10, height = 6, dpi = 300, bg = "white")
print(p1)


# ============================================================================
# 3. PREPARE DATA MATRIX
# ============================================================================

cat("Preparing data matrix for DFA...\n")

# First, remove Johnson
both_ts_long_with_q0 <- both_ts_long_with_q0 %>%
  filter(mgmt_river != "Johnson")

# Convert to wide format
proportion_total_wide <- both_ts_long_with_q0 %>%
  filter(!is.na(total_run_prop)) %>%
  select(mgmt_river, Time_Period, total_run_prop) %>%
  pivot_wider(names_from = Time_Period, values_from = total_run_prop, values_fill = 0)

# Convert to matrix
river_names <- proportion_total_wide$mgmt_river
proportion_total_matrix <- as.matrix(proportion_total_wide[, -1])
rownames(proportion_total_matrix) <- river_names

# Remove rivers with zero production or variance
row_sums <- rowSums(proportion_total_matrix, na.rm = TRUE)
row_vars <- apply(proportion_total_matrix, 1, var, na.rm = TRUE)
good_rivers <- which(row_sums > 0 & row_vars > 1e-10 & !is.na(row_vars))
proportion_total_matrix_clean <- proportion_total_matrix[good_rivers, ]

# Scale data to 0-1 range
proportion_total_matrix_scaled <- t(apply(proportion_total_matrix_clean, 1, function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if(x_max == x_min) return(rep(0.5, length(x)))
  return((x - x_min) / (x_max - x_min))
}))

# Store scaling parameters
scaling_params <- data.frame(
  river = rownames(proportion_total_matrix_clean),
  min_val = apply(proportion_total_matrix_clean, 1, min, na.rm = TRUE),
  max_val = apply(proportion_total_matrix_clean, 1, max, na.rm = TRUE)
)

cat(sprintf("Final data: %d rivers x %d time periods\n", 
            nrow(proportion_total_matrix_scaled), ncol(proportion_total_matrix_scaled)))
cat("Data range after scaling: [", round(min(proportion_total_matrix_scaled, na.rm = TRUE), 3), 
    ", ", round(max(proportion_total_matrix_scaled, na.rm = TRUE), 3), "]\n")

# Create faceted plot of scaled data
scaled_data_long <- proportion_total_matrix_scaled %>%
  as.data.frame() %>%
  mutate(mgmt_river = rownames(.)) %>%
  pivot_longer(cols = -mgmt_river, names_to = "time_period", values_to = "scaled_value") %>%
  mutate(
    year = as.numeric(str_extract(time_period, "\\d{4}")),
    quarter = str_extract(time_period, "Q\\d"),
    time_numeric = case_when(
      quarter == "Q0" ~ year - 0.25,
      quarter == "Q1" ~ year + 0.00,
      quarter == "Q2" ~ year + 0.25,
      quarter == "Q3" ~ year + 0.50,
      quarter == "Q4" ~ year + 0.75
    )
  )

p_scaled <- ggplot(scaled_data_long, aes(x = time_numeric, y = scaled_value)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.5, alpha = 0.7) +
  facet_wrap(~mgmt_river, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = unique(scaled_data_long$year)) +
  labs(
    title = "Scaled Total Run Proportion Data by Management River",
    subtitle = "Data scaled to 0-1 range for each river",
    x = "Year-Quarter",
    y = "Scaled Value (0-1)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIGURE_PATH, "scaled_data_faceted_plot.png"), 
       p_scaled, width = 12, height = 10, dpi = 300, bg = "white")
print(p_scaled)
# ============================================================================
# 4. DFA MODEL COMPARISON Covariates vs no 
# ============================================================================

# Reshape to wide format (columns = Year_Quartile, value = cpue_prop_in_quartile)
cpue_wide <- both_ts_long_with_q0 %>%
  mutate(year_q = paste0(year, "_", quartile)) %>%
  select(year_q, cpue_prop_in_quartile) %>%
  distinct() %>%
  pivot_wider(names_from = year_q, values_from = cpue_prop_in_quartile)

# Convert to 1-row matrix
cpue_cov <- as.matrix(cpue_wide)

runsize_data<- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Escapements_runsize.csv")

runsize_data <- runsize_data[rep(1:nrow(runsize_data), each = 5), ]  # 5 for Q0, Q1, Q2, Q3, Q4

runsize_cov <- matrix(runsize_data$Total.Run, nrow = 1)

all_covariates <- rbind(cpue_cov, runsize_cov)

# Z-score normalize covariates
cpue_cov <- scale(t(cpue_cov))
cpue_cov <- t(cpue_cov)

runsize_cov <- scale(t(runsize_cov))
runsize_cov <- t(runsize_cov)

all_covariates <- rbind(cpue_cov, runsize_cov)

# Fit a 4-state DFA model with all covariates

all_cvs <- MARSS(proportion_total_matrix_scaled,
                 model = list(m = 4, R = "diagonal and unequal"),
                 z.score = TRUE,
                 form = "dfa",
                 control = list(maxit = 50000), 
                 method = "BFGS", 
                 covariates = all_covariates,
                 silent = TRUE)

cpue_only <- MARSS(proportion_total_matrix_scaled,
                   model = list(m = 4, R = "diagonal and unequal"),
                   z.score = TRUE,
                   form = "dfa",
                   control = list(maxit = 50000), 
                   method = "BFGS", 
                   covariates = cpue_cov,
                   silent = TRUE)

runsize_only <- MARSS(proportion_total_matrix_scaled,
                      model = list(m = 4, R = "diagonal and unequal"),
                      z.score = TRUE,
                      form = "dfa",
                      control = list(maxit = 50000), 
                      method = "BFGS", 
                      covariates = runsize_cov,
                      silent = TRUE)

no_covariates <- MARSS(proportion_total_matrix_scaled,
                       model = list(m = 4, R = "diagonal and equal"),
                       z.score = TRUE,
                       form = "dfa",
                       control = list(maxit = 50000), 
                       method = "BFGS",
                       silent = TRUE)

# Extract and compare aic 
aic_values <- c(
  all_cvs = all_cvs$AIC,
  cpue_only = cpue_only$AIC,
  runsize_only = runsize_only$AIC,
  no_covariates = no_covariates$AIC
)

# Print the AIC values
print(aic_values)

# According to AIC, select best model for further analysis

# ============================================================================
# 4. DFA MODEL COMPARISON (1-4 STATES)
# ============================================================================

cat("\n=== COMPARING DFA MODELS (1-4 STATES) ===\n")

# Initialize results storage
model_results <- data.frame(
  states = 1:4,
  rmse = numeric(4),
  aic = numeric(4),
  aicc = numeric(4),
  loglik = numeric(4),
  converged = logical(4)
)

# Test each model
for(m_states in 1:4) {
  cat(sprintf("Fitting %d-state model...\n", m_states))
  
  tryCatch({
    dfa_model <- MARSS(proportion_total_matrix,
                       model = list(m = m_states, R = "diagonal and unequal"),
                       z.score = FALSE,
                       form = "dfa",
                       control = list(maxit = 50000), 
                       covariates = cpue_cov,
                       method = "BFGS", 
                       silent = TRUE)
    
    if(dfa_model$convergence == 0) {
      # Calculate RMSE
      fitted_values <- if(m_states == 1) {
        coef(dfa_model, type = "matrix")$Z %*% dfa_model$states
      } else {
        H_inv <- varimax(coef(dfa_model, type = "matrix")$Z)$rotmat
        coef(dfa_model, type = "matrix")$Z %*% H_inv %*% dfa_model$states
      }
      rmse <- sqrt(mean((proportion_total_matrix - fitted_values)^2, na.rm = TRUE))
      
      # Store results
      model_results[m_states, ] <- list(
        states = m_states,
        rmse = rmse,
        aic = dfa_model$AIC,
        aicc = dfa_model$AICc,
        loglik = dfa_model$logLik,
        converged = TRUE
      )
      
      cat(sprintf("  RMSE: %.4f | AICc: %.2f\n", rmse, dfa_model$AICc))
      
    } else {
      model_results[m_states, "converged"] <- FALSE
      cat("  Model did not converge\n")
    }
    
  }, error = function(e) {
    model_results[m_states, "converged"] <- FALSE
    cat(sprintf("  Error: %s\n", e$message))
  })
}

# Display results
cat("\n=== MODEL COMPARISON RESULTS ===\n")
print(model_results[model_results$converged, ])

# Select best model
best_n_states <- 3
cat(sprintf("\nUsing %d-state model for final analysis\n", best_n_states))

# ============================================================================
# 4. DFA MODEL COMPARISON Different variance structures
# ============================================================================

r_equal <- MARSS(proportion_total_matrix,
                 model = list(m = 2, R = "diagonal and equal"),
                 z.score = TRUE,
                 form = "dfa",
                 control = list(maxit = 50000), 
                 method = "BFGS",
                 silent = TRUE)


r_diagonal <- MARSS(proportion_total_matrix,
                    model = list(m = 2, R = "diagonal and unequal"),
                    z.score = TRUE,
                    form = "dfa",
                    control = list(maxit = 50000), 
                    method = "BFGS",
                    silent = TRUE)


# Extract and compare AIC values
aic_values_variance <- c(
  r_equal = r_equal$AIC,
  r_diagonal = r_diagonal$AIC
)

print(aic_values_variance)

# ============================================================================
# 5. FIT BEST DFA MODEL
## no covariates 
## 2 state model 
## diagonal and unequal variance structure (based on your within quartile results)
# ============================================================================

best_model_total<- MARSS(proportion_total_matrix,
                         model = list(m = 2, R = "diagonal and unequal"),
                         z.score = FALSE,
                         form = "dfa", 
                         covariates = cpue_cov,
                         control = list(maxit = 50000), 
                         method = "BFGS")

# Extract and rotate results
trends_total_orig <- best_model_total$states
Z_total_orig <- coef(best_model_total, type = "matrix")$Z

# Apply varimax rotation
cat("Applying varimax rotation...\n")
varimax_result_total <- varimax(Z_total_orig)
H_inv_total <- varimax_result_total$rotmat

# Rotate loadings and trends
Z_total <- Z_total_orig %*% H_inv_total
trends_total <- solve(H_inv_total) %*% trends_total_orig
colnames(Z_total) <- paste("Trend", 1:2)

# ============================================================================
# 6. CREATE TIME LABELS
# ============================================================================

# Create time labels including Q0
years <- rep(2017:2021, each = 5)  # 5 quarters per year (Q0, Q1, Q2, Q3, Q4)
quarters <- rep(0:4, times = 5)
time_labels <- paste0(years, "-Q", quarters)
time_labels_used <- time_labels[1:ncol(trends_total)]

# ============================================================================
# 7. PLOT DFA RESULTS
# ============================================================================

for (i in 1:2) {
  cat(sprintf("Creating plots for Trend %d...\n", i))
  
  # Trend time series data
  trend_data <- data.frame(
    time = 1:ncol(trends_total),
    value = trends_total[i, ],
    time_label = time_labels_used
  )
  
  # Trend plot
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
  
  # Loadings data
  loading_data <- data.frame(Unit = rownames(Z_total), Loading = Z_total[, i]) %>%
    arrange(desc(abs(Loading)))
  
  # Loadings plot
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
  
  # Combine and save
  combined_plot <- grid.arrange(trend_plot, loadings_plot, ncol = 2, 
                                top = paste("TOTAL RUN PROPORTION - TREND", i))
  
  ggsave(file.path(FIGURE_PATH, paste0("total_run_proportion_trend_", i, "_plot.png")), 
         combined_plot, width = 12, height = 6, dpi = 300, bg = "white")
  
  print(combined_plot)
}

# ============================================================================
# 8. CREATE SPATIAL MAPS
# ============================================================================

cat("\n=== CREATING SPATIAL MAPS ===\n")

# Load spatial data
tryCatch({
  edges <- st_read(SPATIAL_DATA_PATH, quiet = TRUE)
  basin <- st_read(BASIN_PATH, quiet = TRUE)
  
  cat("Spatial data loaded successfully\n")
  
  if ("mgmt_river" %in% colnames(edges)) {
    
    for (trend_num in 1:best_n_states) {
      cat(paste("Creating map for Trend", trend_num, "\n"))
      
      # Trend time series plot
      trend_data <- data.frame(
        time = 1:length(trends_total[trend_num, ]),
        value = as.numeric(trends_total[trend_num, ]),
        time_label = time_labels_used
      )
      
      trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
        geom_line(color = "steelblue", linewidth = 1.5, alpha = 0.8) +
        geom_point(color = "steelblue", size = 3, alpha = 0.9) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
        labs(
          title = paste("Total Run Proportion - Trend", trend_num, "Time Series"),
          x = "Year-Quarter", 
          y = "Trend Value"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
        )
      
      # Prepare spatial data with loadings
      loadings_data <- data.frame(
        mgmt_river = rownames(Z_total),
        loading = Z_total[, trend_num],
        stringsAsFactors = FALSE
      )
      
      edges_with_loadings <- edges %>%
        left_join(loadings_data, by = "mgmt_river") %>%
        filter(!is.na(mgmt_river) & mgmt_river != "") %>%
        mutate(
          stream_order = ifelse(is.na(Str_Order), 3, Str_Order),
          line_width = pmax(0.3, pmin(3.0, 0.3 + (stream_order - min(stream_order, na.rm = TRUE)) * 
                                        (3.0 - 0.3) / (max(stream_order, na.rm = TRUE) - min(stream_order, na.rm = TRUE))))
        )
      
      # Create spatial map
      if (nrow(edges_with_loadings) > 0) {
        if (st_crs(basin) != st_crs(edges_with_loadings)) {
          basin <- st_transform(basin, st_crs(edges_with_loadings))
        }
        
        map_plot <- ggplot() +
          geom_sf(data = basin, fill = "gray95", color = "gray70", 
                  linewidth = 0.5, alpha = 0.3) +
          geom_sf(data = edges_with_loadings, 
                  aes(color = loading, linewidth = stream_order), 
                  alpha = 0.8) +
          scale_color_gradient2(
            low = "blue", mid = "white", high = "red", midpoint = 0,
            name = "Loading\nValue",
            limits = c(-max(abs(loadings_data$loading)), max(abs(loadings_data$loading)))
          ) +
          scale_linewidth_continuous(
            range = c(0.3, 3.0), 
            name = "Stream\nOrder"
          ) +
          coord_sf(datum = NA) +
          labs(
            title = paste("Total Run Proportion - Trend", trend_num, "Loadings Map"),
            subtitle = "Rivers colored by loading values (blue = negative, red = positive)"
          ) +
          theme_void() +
          theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            legend.position = "right",
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)
          )
        
        # Save combined map
        filename <- file.path(MAP_PATH, paste0("total_run_proportion_trend_", trend_num, "_map.png"))
        
        png(filename, width = 16, height = 8, units = "in", res = 300, bg = "white")
        grid.arrange(
          trend_plot, map_plot, 
          ncol = 2, 
          widths = c(1, 1.2),
          top = textGrob(paste("Total Run Proportion - Trend", trend_num, "Analysis"), 
                         gp = gpar(fontsize = 16, fontface = "bold"))
        )
        dev.off()
        
        cat(paste("Saved:", filename, "\n"))
      }
    }
    
  } else {
    cat("Warning: mgmt_river column not found in spatial data\n")
  }
  
}, error = function(e) {
  cat("Error loading spatial data:", e$message, "\n")
})

# ============================================================================
# 9. SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("✓ Data prepared with Q0 entries\n")
cat("✓ DFA models compared (1-4 states)\n")
cat(sprintf("✓ Best model fitted (%d states)\n", best_n_states))
cat("✓ Varimax rotation applied\n")
cat("✓ Trend plots created\n")
cat("✓ Spatial maps generated\n")
cat("\nAll results saved to specified directories.\n")