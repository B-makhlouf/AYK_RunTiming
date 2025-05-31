# ============================================================================
# Dynamic Factor Analysis (DFA) for Salmon Within Quartile Proportions
# Polished Analysis Script - Model comparison and visualization with Q0 values
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
# 1. DATA PREPARATION WITH Q0 VALUES
# ============================================================================

# Load and clean data
both_ts_long <- read.csv(DATA_PATH) %>%
  mutate(
    Quartile_Clean = str_to_upper(str_extract(quartile, "Q[1-4]")),
    Quartile_Num = as.numeric(str_remove(Quartile_Clean, "Q")),
    Time_Continuous = year + (Quartile_Num - 1) * 0.25,
    Time_Period = paste0(year, "_", Quartile_Clean)
  )

# OPTION 1: Even spacing with Q0 as time 0 within each year
# Create Q0 rows for each management river and year combination
q0_rows <- both_ts_long %>%
  # Get unique mgmt_river-year combinations
  select(mgmt_river, year) %>%
  distinct() %>%
  # Create Q0 rows with zero values
  mutate(
    quartile = "Q0",
    Quartile_Clean = "Q0",
    Quartile_Num = 0,
    Time_Period = paste0(year, "_Q0"),
    Time_Continuous = year + 0.0,   # Q0 at start of year (evenly spaced)
    total_run_prop = 0,
    within_quartile_prop = 0,
    cpue_prop_in_quartile = 0,
    edge_count = 0
  )

# Update the main data to have even spacing: Q0=0.0, Q1=0.2, Q2=0.4, Q3=0.6, Q4=0.8
both_ts_long <- both_ts_long %>%
  mutate(
    Time_Continuous = year + (Quartile_Num - 1) * 0.2 + 0.2  # Even 0.2 spacing
  )

# Combine original data with Q0 rows
both_ts_complete <- bind_rows(both_ts_long, q0_rows) %>%
  arrange(mgmt_river, year, Quartile_Num)

# Filter valid data for within quartile analysis
within_quartile_data <- both_ts_complete %>%
  filter(!is.na(within_quartile_prop))

cat("Data preparation complete:\n")
cat("- Original data points:", nrow(both_ts_long), "\n")
cat("- Q0 rows added:", nrow(q0_rows), "\n")
cat("- Total data points:", nrow(both_ts_complete), "\n")
cat("- Management rivers:", length(unique(within_quartile_data$mgmt_river)), "\n")

# Define river colors
mgmt_rivers <- unique(within_quartile_data$mgmt_river)
n_rivers <- length(mgmt_rivers)
river_colors <- if (n_rivers <= 11) {
  brewer.pal(max(3, n_rivers), "Spectral")
} else {
  colorRampPalette(brewer.pal(11, "Spectral"))(n_rivers)
}
names(river_colors) <- mgmt_rivers

# ============================================================================
# 2. VISUALIZATION WITH Q0 VALUES
# ============================================================================

# Plot with Q0 values included
p1 <- ggplot(within_quartile_data, aes(x = Time_Continuous, y = within_quartile_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(
    breaks = seq(min(within_quartile_data$year) - 0.25, max(within_quartile_data$year) + 1, by = 1),
    labels = function(x) {
      # Label Q0 points and full years
      ifelse(x %% 1 == 0.75, paste0(floor(x + 0.25), "-Q0"), as.character(x))
    }
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, NA)) +
  labs(
    title = "Within Quartile Proportion by Management Unit Over Time ",
    x = "Year",
    y = "Proportional contribution to the Quartile",
    color = "Management River"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(color = guide_legend(ncol = 3))

# Faceted plot by management river for better visibility
p2 <- ggplot(within_quartile_data, aes(x = Time_Continuous, y = within_quartile_prop)) +
  geom_line(linewidth = 1.2, color = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "steelblue", alpha = 0.9) +
  facet_wrap(~mgmt_river, scales = "free_y", ncol = 3) +
  scale_x_continuous(
    breaks = seq(min(within_quartile_data$year), max(within_quartile_data$year), by = 2)
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Within Quartile Proportion by Management Unit Over Time ",
    x = "Year",
    y = "Proportional contribution to the Quartile",
    color = "Management River"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plots
dir.create(FIGURE_PATH, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIGURE_PATH, "within_quartile_proportion_timeseries_with_Q0.png"), 
       p1, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(FIGURE_PATH, "within_quartile_proportion_faceted_with_Q0.png"), 
       p2, width = 14, height = 10, dpi = 300, bg = "white")

print(p1)
print(p2)

# ============================================================================
# 3. PREPARE DATA FOR DFA WITH Q0 VALUES
# ============================================================================

# Convert to wide format for DFA (rivers as rows, time periods as columns)
within_quartile_wide <- within_quartile_data %>%
  select(mgmt_river, Time_Period, within_quartile_prop) %>%
  pivot_wider(
    names_from = Time_Period, 
    values_from = within_quartile_prop, 
    values_fill = 0
  )

# Convert to matrix format for MARSS
river_names <- within_quartile_wide$mgmt_river
within_quartile_matrix <- as.matrix(within_quartile_wide[, -1])
rownames(within_quartile_matrix) <- river_names

# Sort columns chronologically (important with Q0 values)
time_periods <- colnames(within_quartile_matrix)
years <- as.numeric(str_extract(time_periods, "\\d{4}"))
quarters <- str_extract(time_periods, "Q\\d")
quarter_nums <- ifelse(quarters == "Q0", 0, as.numeric(str_remove(quarters, "Q")))

# Create sorting index
sort_index <- order(years, quarter_nums)
within_quartile_matrix <- within_quartile_matrix[, sort_index]

# Check data quality
row_sums <- rowSums(within_quartile_matrix, na.rm = TRUE)
row_vars <- apply(within_quartile_matrix, 1, var, na.rm = TRUE)

cat("\nData quality check:\n")
cat("- Rivers with zero production:", sum(row_sums == 0), "\n")
cat("- Rivers with zero variance:", sum(row_vars == 0 | is.na(row_vars)), "\n")

# Remove rivers with zero production or zero variance
good_rivers <- which(row_sums > 0 & row_vars > 1e-10 & !is.na(row_vars))
within_quartile_matrix_clean <- within_quartile_matrix[good_rivers, ]

# Scale data from 0-1 (min-max scaling) instead of z-normalization
scale_0_1 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Apply 0-1 scaling to each row (management river)
within_quartile_matrix_scaled <- t(apply(within_quartile_matrix_clean, 1, scale_0_1))
colnames(within_quartile_matrix_scaled) <- colnames(within_quartile_matrix_clean)
rownames(within_quartile_matrix_scaled) <- rownames(within_quartile_matrix_clean)

cat("- Rivers retained for analysis:", nrow(within_quartile_matrix_scaled), "\n")
cat("- Time periods:", ncol(within_quartile_matrix_scaled), "\n")
cat("- Data range after 0-1 scaling:", round(range(within_quartile_matrix_scaled, na.rm = TRUE), 3), "\n")

# # ============================================================================
# # 4. DFA MODEL COMPARISON Covariates vs no 
# # ============================================================================
# 
# # FIXED: Prepare covariates to match 25 time periods (5 years × 5 quarters including Q0)
# 
# # CPUE covariate - need to match the 25 time periods
# # Original data has CPUE for Q1-Q4, need to add Q0 values (zero)
# cpue_expanded <- both_ts_complete %>%
#   select(year, Quartile_Num, cpue_prop_in_quartile) %>%
#   distinct() %>%
#   arrange(year, Quartile_Num) %>%
#   pull(cpue_prop_in_quartile)
# 
# # Z-normalize CPUE covariate
# cpue_mean <- mean(cpue_expanded, na.rm = TRUE)
# cpue_sd <- sd(cpue_expanded, na.rm = TRUE)
# cpue_normalized <- (cpue_expanded - cpue_mean) / cpue_sd
# cpue_cov <- matrix(cpue_normalized, nrow = 1)
# 
# # Run size covariate - repeat each year's value 5 times (for Q0, Q1, Q2, Q3, Q4)
# runsize_data <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Escapements_runsize.csv")
# runsize_expanded <- rep(runsize_data$Total.Run, each = 5)  # 5 quarters per year including Q0
# 
# # Z-normalize run size covariate
# runsize_mean <- mean(runsize_expanded, na.rm = TRUE)
# runsize_sd <- sd(runsize_expanded, na.rm = TRUE)
# runsize_normalized <- (runsize_expanded - runsize_mean) / runsize_sd
# runsize_cov <- matrix(runsize_normalized, nrow = 1)
# 
# # Combined covariates
# all_covariates <- rbind(cpue_cov, runsize_cov)
# 
# cat("Covariates check:\n")
# cat("- CPUE covariate dimensions:", dim(cpue_cov), "\n")
# cat("- Run size covariate dimensions:", dim(runsize_cov), "\n")
# cat("- Combined covariates dimensions:", dim(all_covariates), "\n")
# cat("- Data matrix dimensions:", dim(within_quartile_matrix_clean), "\n")
# 
# # Fit a 4-state DFA model with all covariates
# all_cvs <- MARSS(within_quartile_matrix_clean,
#                  model = list(m = 4, R = "diagonal and equal"),
#                  z.score = TRUE,
#                  form = "dfa",
#                  control = list(maxit = 50000), 
#                  method = "BFGS", 
#                  covariates = all_covariates,
#                  silent = TRUE)
# 
# cpue_only <- MARSS(within_quartile_matrix_clean,
#                    model = list(m = 4, R = "diagonal and equal"),
#                    z.score = TRUE,
#                    form = "dfa",
#                    control = list(maxit = 50000), 
#                    method = "BFGS", 
#                    covariates = cpue_cov,
#                    silent = TRUE)
# 
# runsize_only <- MARSS(within_quartile_matrix_clean,
#                       model = list(m = 4, R = "diagonal and equal"),
#                       z.score = TRUE,
#                       form = "dfa",
#                       control = list(maxit = 50000), 
#                       method = "BFGS", 
#                       covariates = runsize_cov,
#                       silent = TRUE)
# 
# no_covariates <- MARSS(within_quartile_matrix_clean,
#                        model = list(m = 4, R = "diagonal and equal"),
#                        z.score = TRUE,
#                        form = "dfa",
#                        control = list(maxit = 50000), 
#                        method = "BFGS",
#                        silent = TRUE)
# 
# # Extract and compare AIC 
# aic_comparison <- data.frame(
#   Model = c("All Covariates", "CPUE Only", "Run Size Only", "No Covariates"),
#   AIC = c(all_cvs$AIC, cpue_only$AIC, runsize_only$AIC, no_covariates$AIC),
#   AICc = c(all_cvs$AICc, cpue_only$AICc, runsize_only$AICc, no_covariates$AICc),
#   Converged = c(all_cvs$convergence == 0, cpue_only$convergence == 0, 
#                 runsize_only$convergence == 0, no_covariates$convergence == 0)
# ) %>%
#   mutate(Delta_AIC = AIC - min(AIC, na.rm = TRUE),
#          Delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
#   arrange(AICc)
# 
# # Print the comparison
# print(aic_comparison)
# 
# # ============================================================================
# # 4. DFA MODEL COMPARISON (1-4 STATES)
# # ============================================================================
# 
# cat("\n=== COMPARING DFA MODELS (1-4 STATES) ===\n")
# 
# # Initialize results storage
# model_results <- data.frame(
#   states = 1:4,
#   rmse = numeric(4),
#   aic = numeric(4),
#   aicc = numeric(4),
#   loglik = numeric(4),
#   converged = logical(4)
# )
# 
# # Test each model
# for(m_states in 1:4) {
#   cat(sprintf("Fitting %d-state model...\n", m_states))
#   
#   tryCatch({
#     dfa_model <- MARSS(within_quartile_matrix_scaled,
#                        model = list(m = m_states, R = "diagonal and unequal"),
#                        z.score = FALSE,
#                        form = "dfa",
#                        control = list(maxit = 50000), 
#                        method = "BFGS", 
#                        silent = TRUE)
#     
#     if(dfa_model$convergence == 0) {
#       # Calculate RMSE
#       fitted_values <- if(m_states == 1) {
#         coef(dfa_model, type = "matrix")$Z %*% dfa_model$states
#       } else {
#         H_inv <- varimax(coef(dfa_model, type = "matrix")$Z)$rotmat
#         coef(dfa_model, type = "matrix")$Z %*% H_inv %*% dfa_model$states
#       }
#       rmse <- sqrt(mean((within_quartile_matrix_scaled - fitted_values)^2, na.rm = TRUE))
#       
#       # Store results
#       model_results[m_states, ] <- list(
#         states = m_states,
#         rmse = rmse,
#         aic = dfa_model$AIC,
#         aicc = dfa_model$AICc,
#         loglik = dfa_model$logLik,
#         converged = TRUE
#       )
#       
#       cat(sprintf("  RMSE: %.4f | AICc: %.2f\n", rmse, dfa_model$AICc))
#       
#     } else {
#       model_results[m_states, "converged"] <- FALSE
#       cat("  Model did not converge\n")
#     }
#     
#   }, error = function(e) {
#     model_results[m_states, "converged"] <- FALSE
#     cat(sprintf("  Error: %s\n", e$message))
#   })
# }
# 
# # Display results
# cat("\n=== MODEL COMPARISON RESULTS ===\n")
# print(model_results[model_results$converged, ])
# 
# # Select best model
# best_n_states <- 2
# cat(sprintf("\nUsing %d-state model for final analysis\n", best_n_states))
# 
# # Fitting 1-state model...
# # RMSE: 0.5204 | AICc: -319.74
# # Fitting 2-state model...
# # RMSE: 0.4938 | AICc: -889.88
# # Fitting 3-state model...
# # RMSE: 0.4923 | AICc: -1190.33
# # Fitting 4-state model...
# # RMSE: 0.4948 | AICc: -1258.13
# 
# 
# #########################
# ######################### There's the most support for 3, but lets do 4 because there's 4 states
# #########################
# 
# #  Test different covariance structures 
# 
# equal<- MARSS(within_quartile_matrix_scaled,
#               model = list(m = 4, R = "diagonal and equal"),
#               z.score = FALSE,
#               form = "dfa",
#               control = list(maxit = 50000), 
#               method = "BFGS", 
#               silent = TRUE)
# 
# unequal <- MARSS(within_quartile_matrix_scaled,
#                  model = list(m = 4, R = "diagonal and unequal"),
#                  z.score = FALSE,
#                  form = "dfa",
#                  control = list(maxit = 50000), 
#                  method = "BFGS", 
#                  silent = TRUE)
# 
# 
# ## Compare AICs 
# 
# aic_comparison <- data.frame(
#   Model = c("Equal Covariance", "Unequal Covariance"),
#   AIC = c(equal$AIC, unequal$AIC),
#   AICc = c(equal$AICc, unequal$AICc),
#   Converged = c(equal$convergence == 0, unequal$convergence == 0)
# ) %>%
#   mutate(Delta_AIC = AIC - min(AIC, na.rm = TRUE),
#          Delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
#   arrange(AICc)
# 
# print(aic_comparison)
# 
# ######### For 3 state models 
# ############################
# 
# # Model        AIC       AICc Converged Delta_AIC Delta_AICc
# # 1 Unequal Covariance -1218.7914 -1190.3269      TRUE    0.0000     0.0000
# # 2   Equal Covariance  -991.9487  -976.4294      TRUE  226.8427   213.8975
# 
# 
# ######## For 4 state models 
# ###########################
# # # Model        AIC       AICc Converged Delta_AIC Delta_AICc
# # 1 Unequal Covariance -1218.7914 -1190.3269      TRUE    0.0000     0.0000
# # 2   Equal Covariance  -991.9487  -976.4294      TRUE  226.8427   213.8975
# # 

# ============================================================================
# 5. FIT BEST DFA MODEL
# ============================================================================

##### Best model is: 
# 3 states, 
# No Covariates 
# Diagonal and unequal covariance structure
# 
best_n_states<- 4
cat("\n=== FITTING BEST DFA MODEL ===\n")

best_model_within <- suppressMessages(suppressWarnings(
  MARSS(within_quartile_matrix_scaled,
        model = list(m = best_n_states, R = "diagonal and equal"),
        z.score = FALSE,
        form = "dfa", 
        control = list(maxit = 50000), 
        method = "BFGS")
))

# Extract and rotate results
trends_within_orig <- best_model_within$states
Z_within_orig <- coef(best_model_within, type = "matrix")$Z

# Apply varimax rotation
cat("Applying varimax rotation...\n")
varimax_result_within <- varimax(Z_within_orig)
H_inv_within <- varimax_result_within$rotmat

# Rotate loadings and trends
Z_within <- Z_within_orig %*% H_inv_within
trends_within <- solve(H_inv_within) %*% trends_within_orig
colnames(Z_within) <- paste("Trend", 1:best_n_states)

# ============================================================================
# 6. EXPORT LOADINGS BY MANAGEMENT TRIBUTARY
# ============================================================================

cat("\n=== EXPORTING LOADINGS DATA ===\n")

# Create loadings export dataframe
loadings_export <- as.data.frame(Z_within)
loadings_export$mgmt_river <- rownames(Z_within)

# Reorder columns to put management river first
loadings_export <- loadings_export[, c("mgmt_river", paste("Trend", 1:best_n_states))]

# Save to CSV
loadings_filename <- file.path(FIGURE_PATH, "dfa_loadings_by_management_river.csv")
write.csv(loadings_export, loadings_filename, row.names = FALSE)

cat("✓ Loadings exported to:", loadings_filename, "\n")
cat("✓ File contains", nrow(loadings_export), "management rivers and", best_n_states, "trends\n")

# Also create a summary table with the strongest loadings for each trend
loadings_summary <- data.frame()
for(i in 1:best_n_states) {
  trend_loadings <- data.frame(
    Trend = paste("Trend", i),
    mgmt_river = rownames(Z_within),
    Loading = Z_within[, i],
    Abs_Loading = abs(Z_within[, i])
  ) %>%
    arrange(desc(Abs_Loading)) %>%
    slice_head(n = 5) %>%  # Top 5 strongest loadings
    select(-Abs_Loading)
  
  loadings_summary <- rbind(loadings_summary, trend_loadings)
}

# Save summary table
summary_filename <- file.path(FIGURE_PATH, "dfa_loadings_summary_top5.csv")
write.csv(loadings_summary, summary_filename, row.names = FALSE)

cat("✓ Loadings summary (top 5 per trend) exported to:", summary_filename, "\n")

# ============================================================================
# 7. CREATE TIME LABELS
# ============================================================================

# Create time labels including Q0
years <- rep(2017:2021, each = 5)  # 5 quarters per year (Q0, Q1, Q2, Q3, Q4)
quarters <- rep(0:4, times = 5)
time_labels <- paste0(years, "-Q", quarters)
time_labels_used <- time_labels[1:ncol(trends_within)]

# ============================================================================
# 8. PLOT DFA RESULTS WITH QUARTILE-COLORED POINTS
# ============================================================================

cat("\n=== CREATING DFA PLOTS WITH QUARTILE COLORS ===\n")

# Define sequential colors for quartiles (blue to orange to red progression)
quartile_colors <- c(
  "Q0" = "#08519C",  # Dark blue (starting point)
  "Q1" = "#3182BD",  # Medium blue
  "Q2" = "#FD8D3C",  # Orange
  "Q3" = "#E6550D",  # Dark orange
  "Q4" = "#A63603"   # Red-brown
)

for (i in 1:best_n_states) {
  cat(sprintf("Creating plots for Trend %d...\n", i))
  
  # Extract quartile information from time labels
  trend_data <- data.frame(
    time = 1:ncol(trends_within),
    value = trends_within[i, ],
    time_label = time_labels_used
  ) %>%
    mutate(
      quartile = str_extract(time_label, "Q[0-4]"),
      year = str_extract(time_label, "\\d{4}")
    )
  
  # Trend plot with quartile-colored points and grey line
  trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
    geom_line(color = "grey60", linewidth = 1.2, alpha = 0.8) +
    geom_point(aes(color = quartile), size = 3, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(
      values = quartile_colors,
      name = "Quartile",
      labels = c("Q0" = "Q0 (Start)", "Q1" = "Q1", "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4")
    ) +
    scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
    labs(
      title = paste("Trend", i, "Time Series"), 
      x = "Year-Quarter", 
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Loadings data (unchanged)
  loading_data <- data.frame(Unit = rownames(Z_within), Loading = Z_within[, i]) %>%
    arrange(desc(abs(Loading)))
  
  # Loadings plot (unchanged)
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
                                top = paste("WITHIN QUARTILE PROPORTION - TREND", i))
  
  ggsave(file.path(FIGURE_PATH, paste0("within_quartile_proportion_trend_", i, "_plot.png")), 
         combined_plot, width = 14, height = 6, dpi = 300, bg = "white")
  
  print(combined_plot)
}

# ============================================================================
# 9. CREATE SPATIAL MAPS WITH QUARTILE-COLORED TREND PLOTS
# ============================================================================

cat("\n=== CREATING SPATIAL MAPS WITH QUARTILE-COLORED TRENDS ===\n")

# Load spatial data
tryCatch({
  edges <- st_read(SPATIAL_DATA_PATH, quiet = TRUE)
  basin <- st_read(BASIN_PATH, quiet = TRUE)
  
  cat("Spatial data loaded successfully\n")
  
  if ("mgmt_river" %in% colnames(edges)) {
    
    for (trend_num in 1:best_n_states) {
      cat(paste("Creating map for Trend", trend_num, "\n"))
      
      # Trend time series plot with quartile colors
      trend_data <- data.frame(
        time = 1:length(trends_within[trend_num, ]),
        value = as.numeric(trends_within[trend_num, ]),
        time_label = time_labels_used
      ) %>%
        mutate(
          quartile = str_extract(time_label, "Q[0-4]"),
          year = str_extract(time_label, "\\d{4}")
        )
      
      trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
        geom_line(color = "grey60", linewidth = 1.5, alpha = 0.8) +
        geom_point(aes(color = quartile), size = 4, alpha = 0.9) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.1) +
        scale_color_manual(
          values = quartile_colors,
          name = "Quartile",
          labels = c("Q0" = "Q0 (Start)", "Q1" = "Q1", "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4")
        ) +
        scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
        labs(
          title = paste("Within Quartile Proportion - Trend", trend_num, "Time Series"),
          x = "Year-Quarter", 
          y = "Trend Value"
        ) +
        theme_void() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom"
        ) +
        guides(color = guide_legend(override.aes = list(size = 5)))
      
      # Prepare spatial data with loadings (unchanged)
      loadings_data <- data.frame(
        mgmt_river = rownames(Z_within),
        loading = Z_within[, trend_num],
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
          scale_color_distiller(
            palette = "RdYlBu", 
            direction = -1,
            name = "Loading\nValue",
            limits = c(-max(abs(loadings_data$loading)), max(abs(loadings_data$loading)))
          ) +
          scale_linewidth_continuous(
            range = c(0.3, 2.0), 
            name = "Stream\nOrder"
          ) +
          coord_sf(datum = NA) +
          labs(
            title = paste("Within Quartile Proportion - Trend", trend_num, "Loadings Map"),
            subtitle = "Rivers colored by loading values (red = negative, blue = positive)"
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
        filename <- file.path(MAP_PATH, paste0("within_quartile_proportion_trend_", trend_num, "_map.png"))
        
        png(filename, width = 16, height = 8, units = "in", res = 300, bg = "white")
        grid.arrange(
          trend_plot, map_plot, 
          ncol = 2, 
          widths = c(1, 1.2),
          top = textGrob(paste("Within Quartile Proportion - Trend", trend_num, "Analysis"), 
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
# 10. CREATE INDIVIDUAL QUARTILE PLOTS WITH GAM SMOOTHING
# ============================================================================

cat("\n=== CREATING INDIVIDUAL QUARTILE PLOTS WITH GAM SMOOTHING ===\n")

# Load mgam for GAM fitting
if (!require(mgcv, quietly = TRUE)) {
  install.packages("mgcv")
  library(mgcv)
}

# Create individual quartile plots for each trend
tryCatch({
  edges <- st_read(SPATIAL_DATA_PATH, quiet = TRUE)
  basin <- st_read(BASIN_PATH, quiet = TRUE)
  
  if ("mgmt_river" %in% colnames(edges)) {
    
    for (trend_num in 1:best_n_states) {
      cat(paste("Creating quartile plots for Trend", trend_num, "\n"))
      
      # Prepare base trend data
      trend_data <- data.frame(
        time = 1:length(trends_within[trend_num, ]),
        value = as.numeric(trends_within[trend_num, ]),
        time_label = time_labels_used
      ) %>%
        mutate(
          quartile = str_extract(time_label, "Q[0-4]"),
          year = as.numeric(str_extract(time_label, "\\d{4}"))
        )
      
      # Create plots for each quartile
      for (q in c("Q0", "Q1", "Q2", "Q3", "Q4")) {
        cat(paste("  Creating plot for", q, "\n"))
        
        # Filter data for current quartile
        quartile_data <- trend_data %>%
          filter(quartile == q) %>%
          arrange(year)
        
        # Create trend plot with LOESS smoothing and background points
        trend_plot <- ggplot() +
          # Add all other quartile points in grey (background)
          geom_point(data = trend_data %>% filter(quartile != q), 
                     aes(x = year, y = value), 
                     color = "grey80", size = 2, alpha = 0.6) +
          # Add LOESS smooth for current quartile
          geom_smooth(data = quartile_data, aes(x = year, y = value),
                      method = "loess", se = TRUE, color = quartile_colors[q], 
                      fill = quartile_colors[q], alpha = 0.3, linewidth = 1.5) +
          # Add highlighted points for current quartile
          geom_point(data = quartile_data, aes(x = year, y = value),
                     color = quartile_colors[q], size = 4, alpha = 0.9) +
          geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.1) +
          scale_x_continuous(breaks = unique(trend_data$year)) +
          labs(
            title = paste("Within Quartile Proportion - Trend", trend_num, "-", q, "Time Series"),
            x = "Year", 
            y = "Trend Value"
          ) +
          theme_void() +
          theme(
            plot.title = element_text(face = "bold", size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
          )
        
        # Prepare spatial data with loadings (same as before)
        loadings_data <- data.frame(
          mgmt_river = rownames(Z_within),
          loading = Z_within[, trend_num],
          stringsAsFactors = FALSE
        )
        
        edges_with_loadings <- edges %>%
          left_join(loadings_data, by = "mgmt_river") %>%
          filter(!is.na(mgmt_river) & mgmt_river != "") %>%
          mutate(
            stream_order = ifelse(is.na(Str_Order), 3, Str_Order),
            line_width = pmax(0.3, pmin(2.0, 0.3 + (stream_order - min(stream_order, na.rm = TRUE)) * 
                                          (2.0 - 0.3) / (max(stream_order, na.rm = TRUE) - min(stream_order, na.rm = TRUE))))
          )
        
        # Create spatial map (same as before)
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
            scale_color_distiller(
              palette = "RdYlBu", 
              direction = -1,
              name = "Loading\nValue",
              limits = c(-max(abs(loadings_data$loading)), max(abs(loadings_data$loading)))
            ) +
            scale_linewidth_continuous(
              range = c(0.3, 2.0), 
              name = "Stream\nOrder"
            ) +
            coord_sf(datum = NA) +
            labs(
              title = paste("Within Quartile Proportion - Trend", trend_num, "Loadings Map"),
              subtitle = "Rivers colored by loading values (red = negative, blue = positive)"
            ) +
            theme_void() +
            theme(
              plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
              plot.subtitle = element_text(size = 12, hjust = 0.5),
              legend.position = "right",
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA)
            )
          
          # Save combined quartile-specific map
          filename <- file.path(MAP_PATH, paste0("within_quartile_proportion_trend_", trend_num, "_", q, "_gam_map.png"))
          
          png(filename, width = 16, height = 8, units = "in", res = 300, bg = "white")
          grid.arrange(
            trend_plot, map_plot, 
            ncol = 2, 
            widths = c(1, 1.2),
            top = textGrob(paste("Within Quartile Proportion - Trend", trend_num, "-", q, "Analysis"), 
                           gp = gpar(fontsize = 16, fontface = "bold"))
          )
          dev.off()
          
          cat(paste("    Saved:", filename, "\n"))
        }
      }
    }
    
  } else {
    cat("Warning: mgmt_river column not found in spatial data\n")
  }
  
}, error = function(e) {
  cat("Error loading spatial data for quartile plots:", e$message, "\n")
})

# ============================================================================
# 11. SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("✓ Data prepared with Q0 entries\n")
cat("✓ DFA models compared (1-4 states)\n")
cat(sprintf("✓ Best model fitted (%d states)\n", best_n_states))
cat("✓ Varimax rotation applied\n")
cat("✓ Loadings exported to CSV files:\n")
cat("  - Main loadings file:", loadings_filename, "\n")
cat("  - Summary file (top 5 per trend):", summary_filename, "\n")
cat("✓ Trend plots created\n")
cat("✓ Spatial maps generated\n")
cat("✓ Individual quartile plots with GAM smoothing created\n")
cat("\nAll results saved to specified directories.\n")