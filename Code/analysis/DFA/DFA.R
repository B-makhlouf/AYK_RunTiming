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
  library(scales)
})

# Set paths
DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
BASE_FIGURE_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Figures"
BASE_MAP_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/DFA/Maps"
SPATIAL_DATA_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"

# Create organized directory structure
FIGURE_PATH <- file.path(BASE_FIGURE_PATH, "total_run_prop")
MAP_PATH <- file.path(BASE_MAP_PATH, "total_run_prop")
dir.create(FIGURE_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(MAP_PATH, recursive = TRUE, showWarnings = FALSE)

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

cat(sprintf("Data loaded: %d rows\n", nrow(both_ts_long)))

# Remove Johnson river and filter valid data
both_ts_long <- both_ts_long %>%
  filter(mgmt_river != "Johnson")

proportion_total_data <- both_ts_long %>%
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

# Plot initial time series
p1 <- ggplot(proportion_total_data, aes(x = Time_Continuous, y = total_run_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = river_colors) +
  scale_x_continuous(breaks = unique(proportion_total_data$year)) +
  scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, NA)) +
  labs(
    title = "Proportion of Total by Management River Over Time",
    x = "Year (Quarterly)",
    y = "Proportion of Total (%)",
    color = "Management River"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(FIGURE_PATH, "total_run_prop_timeseries.png"), 
       p1, width = 12, height = 8, dpi = 300, bg = "white")
print(p1)

# ============================================================================
# 2. PREPARE DATA MATRIX FOR DFA
# ============================================================================
cat("Preparing data matrix for DFA...\n")

# Convert to wide format
proportion_total_wide <- both_ts_long %>%
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

cat(sprintf("Matrix dimensions: %d rivers × %d time points\n", 
            nrow(proportion_total_matrix_clean), ncol(proportion_total_matrix_clean)))

# Create faceted plot of raw data
raw_data_long <- proportion_total_matrix_clean %>%
  as.data.frame() %>%
  mutate(mgmt_river = rownames(.)) %>%
  pivot_longer(cols = -mgmt_river, names_to = "time_period", values_to = "proportion") %>%
  mutate(
    year = as.numeric(str_extract(time_period, "\\d{4}")),
    quarter = str_extract(time_period, "Q\\d"),
    time_numeric = case_when(
      quarter == "Q1" ~ year + 0.00,
      quarter == "Q2" ~ year + 0.25,
      quarter == "Q3" ~ year + 0.50,
      quarter == "Q4" ~ year + 0.75
    )
  )

p_raw <- ggplot(raw_data_long, aes(x = time_numeric, y = proportion)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.5, alpha = 0.7) +
  facet_wrap(~mgmt_river, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = unique(raw_data_long$year)) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(
    title = "Raw Total Run Proportion Data by Management River",
    x = "Year-Quarter",
    y = "Proportion of Total Run (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(FIGURE_PATH, "raw_data_faceted_plot.png"), 
       p_raw, width = 12, height = 10, dpi = 300, bg = "white")
print(p_raw)

# ============================================================================
# 3. MODEL SELECTION: FIND OPTIMAL NUMBER OF TRENDS
# ============================================================================
cat("Running model selection for optimal number of trends...\n")

# Fit DFA models with different numbers of states
oneState <- MARSS(proportion_total_matrix_scaled,
                  model = list(m = 1, R = "diagonal and unequal"),
                  z.score = FALSE,
                  form = "dfa",
                  control = list(maxit = 50000), 
                  method = "BFGS", 
                  silent = TRUE)

twoState <- MARSS(proportion_total_matrix_scaled,
                  model = list(m = 2, R = "diagonal and unequal"),
                  z.score = FALSE,
                  form = "dfa",
                  control = list(maxit = 50000), 
                  method = "BFGS",
                  silent = TRUE)

threeState <- MARSS(proportion_total_matrix_scaled,
                    model = list(m = 3, R = "diagonal and unequal"),
                    z.score = FALSE,
                    form = "dfa",
                    control = list(maxit = 50000), 
                    method = "BFGS",
                    silent = TRUE)

fourState <- MARSS(proportion_total_matrix_scaled,
                   model = list(m = 4, R = "diagonal and unequal"),
                   z.score = FALSE,
                   form = "dfa",
                   control = list(maxit = 50000), 
                   method = "BFGS",
                   silent = TRUE)

# First, calculate R² for 1-state model
Z1 <- coef(oneState, type = "matrix")$Z
y_hat1 <- Z1 %*% oneState$states
y_obs <- proportion_total_matrix_scaled
SSR1 <- sum((y_obs - y_hat1)^2, na.rm = TRUE)
SST <- sum((y_obs - rowMeans(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
R2_1 <- 1 - (SSR1 / SST)

# Two State Model
Z2 <- coef(twoState, type = "matrix")$Z
y_hat2 <- Z2 %*% twoState$states
SSR2 <- sum((y_obs - y_hat2)^2, na.rm = TRUE)
SST <- sum((y_obs - rowMeans(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
R2_2 <- 1 - (SSR2 / SST)

# Three State Model  
Z3 <- coef(threeState, type = "matrix")$Z
y_hat3 <- Z3 %*% threeState$states
SSR3 <- sum((y_obs - y_hat3)^2, na.rm = TRUE)
R2_3 <- 1 - (SSR3 / SST)

# Four State Model
Z4 <- coef(fourState, type = "matrix")$Z
y_hat4 <- Z4 %*% fourState$states
SSR4 <- sum((y_obs - y_hat4)^2, na.rm = TRUE)
R2_4 <- 1 - (SSR4 / SST)

# Print R² values
cat("Two State R²:", round(R2_2, 3), "\n")
cat("Three State R²:", round(R2_3, 3), "\n") 
cat("Four State R²:", round(R2_4, 3), "\n")

# Calculate Partial R² values
# Partial R² = (R²_full - R²_reduced) / (1 - R²_reduced)

partial_R2_2vs1 <- (R2_2 - R2_1) / (1 - R2_1)  # 2nd state vs 1st state
partial_R2_3vs2 <- (R2_3 - R2_2) / (1 - R2_2)  # 3rd state vs 2nd state  
partial_R2_4vs3 <- (R2_4 - R2_3) / (1 - R2_3)  # 4th state vs 3rd state

model_comparison <- data.frame(
  Model = c("1 State", "2 State", "3 State", "4 State"),
  R_squared = round(c(R2_1, R2_2, R2_3, R2_4), 4),
  Partial_R_squared = c(NA, 
                        round(partial_R2_2vs1, 4),
                        round(partial_R2_3vs2, 4), 
                        round(partial_R2_4vs3, 4)),
  Variance_Explained = round(c(R2_1, R2_2, R2_3, R2_4) * 100, 1),
  Additional_Variance = c(round(R2_1 * 100, 1),
                          round((R2_2 - R2_1) * 100, 1),
                          round((R2_3 - R2_2) * 100, 1),
                          round((R2_4 - R2_3) * 100, 1))
)

# Print results
cat("=== DFA Model Comparison ===\n")
print(model_comparison)

#### BEST MODEL IS 2 states

# ============================================================================
# 4. FIT BEST DFA MODEL WITH VARIMAX ROTATION
# ============================================================================
cat("Fitting best DFA model with varimax rotation...\n")

best_m<- 2
best_model_total <- MARSS(proportion_total_matrix,
                          model = list(m = 2, R = "diagonal and unequal"),
                          z.score = FALSE,
                          form = "dfa", 
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

# Create time labels including Q0
years <- rep(2017:2021, each = 4)  # 4 quarters per year (Q1, Q2, Q3, Q4)
quarters <- rep(1:4, times = 5)
time_labels <- paste0(years, "-Q", quarters)
time_labels_used <- time_labels[1:ncol(trends_total)]


# ============================================================================
# 5. CREATE TREND PLOTS WITH QUARTILE-COLORED POINTS
# ============================================================================
cat("Creating trend visualizations...\n")

# Create color palette for quartiles
quartile_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00")  # Q1=Red, Q2=Blue, Q3=Green, Q4=Orange

for (i in 1:2) {  # Using 2-trend model as determined by model selection
  cat(sprintf("Creating plots for Trend %d...\n", i))
  
  # Trend time series data with quartile information
  trend_data <- data.frame(
    time = 1:ncol(trends_total),
    value = trends_total[i, ],
    time_label = time_labels_used,
    quartile = rep(c("Q1", "Q2", "Q3", "Q4"), length.out = ncol(trends_total))
  )
  
  # Trend plot with quartile-colored points
  trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
    geom_line(color = "black", linewidth = 1.5) +
    geom_point(aes(color = quartile), size = 3, alpha = 0.8) +
    scale_color_manual(
      values = setNames(quartile_colors, c("Q1", "Q2", "Q3", "Q4")),
      name = "Quartile"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6, color = "gray50") +
    scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
    labs(
      title = paste("Trend", i, "Time Series"),
      x = "Year-Quarter", 
      y = "Trend Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      legend.position = "bottom"
    )
  
  # Loadings plot
  loading_data <- data.frame(Unit = rownames(Z_total), Loading = Z_total[, i]) %>%
    arrange(desc(abs(Loading)))
  
  loadings_plot <- ggplot(loading_data, aes(x = reorder(Unit, Loading), 
                                            y = Loading, fill = Loading > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "coral"),
      name = "Loading Direction", 
      labels = c("Negative", "Positive")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6, color = "gray50") +
    labs(
      title = paste("Trend", i, "Loadings"), 
      x = "Management River", 
      y = "Loading Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5), 
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  # Combine and save
  combined_plot <- grid.arrange(
    trend_plot, loadings_plot, 
    ncol = 2,
    top = textGrob(paste("Total Run Proportion - Trend", i), 
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )
  
  ggsave(file.path(FIGURE_PATH, paste0("total_run_proportion_trend_", i, "_plot.png")), 
         combined_plot, width = 14, height = 8, dpi = 300, bg = "white")
  
  print(combined_plot)
}

# ============================================================================
# 6. CREATE SPATIAL MAPS (CORRECTED VERSION)
# ============================================================================
cat("Creating spatial maps...\n")

# Load spatial data
tryCatch({
  edges <- st_read(SPATIAL_DATA_PATH, quiet = TRUE)
  basin <- st_read(BASIN_PATH, quiet = TRUE)
  
  cat("Spatial data loaded successfully\n")
  
  if ("mgmt_river" %in% colnames(edges)) {
    
    for (trend_num in 1:best_m) {
      cat(paste("Creating comprehensive map for Trend", trend_num, "\n"))
      
      # Trend time series plot for spatial map (CORRECTED: using trends_total)
      trend_data <- data.frame(
        time = 1:length(trends_total[trend_num, ]),
        value = as.numeric(trends_total[trend_num, ]),
        time_label = time_labels_used,
        quartile = rep(c("Q1", "Q2", "Q3", "Q4"), length.out = length(trends_total[trend_num, ]))
      )
      
      trend_plot <- ggplot(trend_data, aes(x = time, y = value)) +
        geom_line(color = "black", linewidth = 1.8, alpha = 0.9) +
        geom_point(aes(color = quartile), size = 3.5, alpha = 0.9) +
        scale_color_manual(
          values = setNames(quartile_colors, c("Q1", "Q2", "Q3", "Q4")),
          name = "Quartile"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6, color = "gray50") +
        scale_x_continuous(breaks = trend_data$time, labels = trend_data$time_label) +
        labs(
          title = paste("Trend", trend_num, "Time Series"),
          x = "Year-Quarter", 
          y = "Trend Value"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
          legend.position = "bottom"
        )
      
      # Prepare spatial data with loadings (CORRECTED: using Z_total)
      loadings_data <- data.frame(
        mgmt_river = rownames(Z_total),
        loading = Z_total[, trend_num],
        stringsAsFactors = FALSE
      )
      
      # Join with spatial data
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
        
        # Ensure consistent CRS
        if (st_crs(basin) != st_crs(edges_with_loadings)) {
          basin <- st_transform(basin, st_crs(edges_with_loadings))
        }
        
        # Create the map
        map_plot <- ggplot() +
          geom_sf(data = basin, fill = "gray95", color = "gray70", 
                  linewidth = 0.5, alpha = 0.3) +
          geom_sf(data = edges_with_loadings, 
                  aes(color = loading, linewidth = stream_order), 
                  alpha = 0.8) +
          scale_color_gradient2(
            low = "blue", mid = "white", high = "red", midpoint = 0,
            name = "Loading\nValue",
            limits = c(-max(abs(loadings_data$loading)), max(abs(loadings_data$loading))),
            labels = function(x) sprintf("%.2f", x)
          ) +
          scale_linewidth_continuous(
            range = c(0.3, 3.0), 
            name = "Stream\nOrder"
          ) +
          coord_sf(datum = NA) +
          labs(
            title = paste("Trend", trend_num, "Spatial Loadings"),
            subtitle = "Rivers colored by loading values (blue = negative, red = positive)",
            caption = "Line thickness represents stream order"
          ) +
          theme_void() +
          theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            plot.caption = element_text(size = 10, hjust = 0.5),
            legend.position = "right",
            legend.title = element_text(face = "bold", size = 10),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)
          )
        
        # Save combined figure
        filename <- file.path(MAP_PATH, paste0("total_run_proportion_trend_", trend_num, "_comprehensive_map.png"))
        
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
        
      } else {
        cat(paste("No spatial data available for Trend", trend_num, "\n"))
      }
    }
    
  } else {
    cat("Warning: mgmt_river column not found in spatial data\n")
  }
  
}, error = function(e) {
  cat("Error loading spatial data:", e$message, "\n")
})
