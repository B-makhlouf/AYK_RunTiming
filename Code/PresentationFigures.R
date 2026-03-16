################################################################################
# PRESENTATION FIGURES - KUSKOKWIM WATERSHED CLUSTER ANALYSIS
# Exports: Figures/PresFigures (subdirectory of OUTPUT_DIR)
#
# Contents:
#   PRES-A: Individual-year cumulative distribution plots WITHOUT CPUE line
#   PRES-B: Individual-year cumulative distribution plots WITH CPUE line
#   PRES-C: CPUE foregone vs. protection plot with FLIPPED axes
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(RColorBrewer)
library(scales)
library(cowplot)
library(patchwork)
library(tidyr)
library(grid)

################################################################################
# SETUP AND CONFIGURATION  (unchanged from main script)
################################################################################

OUTPUT_DIR    <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming"
PRES_DIR      <- file.path(OUTPUT_DIR, "Figures", "PresFigures")
CSV_DIR       <- file.path(OUTPUT_DIR, "CSV_Exports")

SPATIAL_PATH      <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH        <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
MGMT_DATA_PATH    <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
CUMULATIVE_DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming/Analysis_Results/Cumulative_Distribution/CSV/cumulative_distribution_data_Kusko_3day_intervals.csv"
CPUE_DATA_PATH    <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming/Analysis_Results/RemainingCPUE/daily_cpue_with_remaining_proportion.csv"

CLOSURE_START_DOY <- 153  # June 1
CLOSURE_END_DOY   <- 162  # June 11

dir.create(PRES_DIR, recursive = TRUE, showWarnings = FALSE)

# ----- Cluster assignments -----
cluster_assignments <- tibble::tribble(
  ~mgmt_river,               ~cluster,
  "N. Fork Kusko",           "Cluster 1",
  "E. Fork Kuskokwim River", "Cluster 1",
  "E. Fork Kuskokwim",       "Cluster 1",
  "S. Fork Kusko",           "Cluster 2",
  "Upper Kusko Main",        "Cluster 2",
  "Big River",               "Cluster 2",
  "Swift",                   "Cluster 2",
  "Tatlawiksuk",             "Cluster 2",
  "Stony",                   "Cluster 3",
  "Takotna and Nixon Fork",  "Cluster 4",
  "Holitna",                 "Cluster 4",
  "George",                  "Cluster 4",
  "Hoholitna",               "Cluster 4",
  "Middle Kusko Main",       "Cluster 5",
  "Oskakawlik",              "Cluster 5",
  "Holokuk",                 "Cluster 5",
  "Kisaralik",               "Cluster 5",
  "Tuluksak",                "Cluster 5",
  "Aniak",                   "Cluster 6",
  "Kwethluk",                "Cluster 6"
)

cluster_order <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6")

cluster_colors <- colorRampPalette(c("#8B0000","#FF0000","#FF8C00","#87CEEB","#1E90FF","#000080"))(6)
names(cluster_colors) <- cluster_order

################################################################################
# DATA LOADING AND PREPARATION  (same logic as main script)
################################################################################

MIN_DOY <- 149
MAX_DOY <- 205

# ----- Cumulative distribution data -----
cumulative_data <- read_csv(CUMULATIVE_DATA_PATH) %>%
  filter(mgmt_river != "Johnson", mgmt_river != "Lower Kusko") %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster)) %>%
  filter(doy >= MIN_DOY & doy <= MAX_DOY)

available_years <- sort(unique(cumulative_data$year))
year_range      <- paste0(min(available_years), "-", max(available_years))

# ----- CPUE data -----
cpue_data_full <- read_csv(CPUE_DATA_PATH)

cpue_data <- cpue_data_full %>%
  filter(year %in% available_years) %>%
  filter(doy >= MIN_DOY & doy <= MAX_DOY) %>%
  mutate(remaining_proportion_pct = remaining_proportion * 100)

# ----- Cluster-level cumulative aggregation (interval-based, corrected) -----
cluster_annual_totals <- cumulative_data %>%
  group_by(year, cluster, mgmt_river) %>%
  summarise(final_production = max(cumulative_production, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, cluster) %>%
  summarise(cluster_annual_total = sum(final_production, na.rm = TRUE), .groups = "drop")

cluster_aggregated <- cumulative_data %>%
  group_by(year, mgmt_river, cluster) %>%
  arrange(doy) %>%
  mutate(interval_production = cumulative_production - lag(cumulative_production, default = 0)) %>%
  ungroup() %>%
  group_by(year, doy, cluster) %>%
  summarise(cluster_interval_production = sum(interval_production, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, cluster) %>%
  arrange(doy) %>%
  mutate(cluster_cumulative_production = cumsum(cluster_interval_production)) %>%
  ungroup() %>%
  left_join(cluster_annual_totals, by = c("year", "cluster")) %>%
  mutate(
    cluster_cumulative_percent = (cluster_cumulative_production / cluster_annual_total) * 100,
    cluster = factor(cluster, levels = cluster_order)
  )

################################################################################
# SHARED PLOT ELEMENTS  (used by both PRES-A and PRES-B)
################################################################################

x_breaks <- seq(152, 205, by = 10)
x_labels <- function(x) format(as.Date(x - 1, origin = "2020-01-01"), "%b %d")

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position  = "none",
    panel.spacing.y  = unit(0.5, "lines"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 14, color = "grey20", face = "bold"),
    axis.text.y      = element_text(size = 14, color = "grey20", face = "bold"),
    axis.title.x     = element_text(size = 16, color = "grey20", face = "bold"),
    axis.title.y     = element_text(size = 16, color = "grey20", face = "bold")
  )

################################################################################
# PRES-A: INDIVIDUAL YEAR PLOTS — WITHOUT CPUE LINE
################################################################################

cat("\n=== PRES-A: Individual cumulative distribution plots (no CPUE) ===\n")

for (yr in available_years) {
  yr_data <- cluster_aggregated %>% filter(year == yr)

  p <- ggplot() +
    geom_hline(yintercept = c(25, 50, 75), color = "gray70",
               linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 162, color = "grey50",
               linetype = "dashed", linewidth = 1.2, alpha = 0.8) +
    geom_line(data  = yr_data,
              aes(x = doy, y = cluster_cumulative_percent, color = cluster),
              linewidth = 1.2, alpha = 0.9) +
    geom_point(data = yr_data,
               aes(x = doy, y = cluster_cumulative_percent, color = cluster),
               size = 0.8, alpha = 0.7) +
    scale_color_manual(
      values = cluster_colors,
      breaks = cluster_order,
      labels = function(x) gsub("Cluster", "Group", x),
      guide  = guide_legend(override.aes = list(linewidth = 2))
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 0), "%"),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(MIN_DOY, MAX_DOY),
      expand = c(0, 0)
    ) +
    labs(
      title = paste("Cumulative distribution vs. Time —", yr),
      x     = "Date",
      y     = "Cumulative Percent of Group's Total Production",
      color = "Group"
    ) +
    base_theme

  fname <- file.path(PRES_DIR, paste0("cumulative_dist_no_cpue_", yr, ".png"))
  ggsave(fname, p, width = 7, height = 5, dpi = 300, bg = "white")
  cat("  Saved:", basename(fname), "\n")
}

cat("=== PRES-A COMPLETE: ", length(available_years), "plots saved ===\n\n")

################################################################################
# PRES-B: INDIVIDUAL YEAR PLOTS — WITH CPUE LINE
################################################################################

cat("=== PRES-B: Individual cumulative distribution plots (with CPUE) ===\n")

for (yr in available_years) {
  yr_data   <- cluster_aggregated %>% filter(year == yr)
  yr_cpue   <- cpue_data          %>% filter(year == yr)

  p <- ggplot() +
    geom_hline(yintercept = c(25, 50, 75), color = "gray70",
               linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 162, color = "grey50",
               linetype = "dashed", linewidth = 1.2, alpha = 0.8) +
    geom_line(data  = yr_data,
              aes(x = doy, y = cluster_cumulative_percent, color = cluster),
              linewidth = 1.2, alpha = 0.9) +
    geom_point(data = yr_data,
               aes(x = doy, y = cluster_cumulative_percent, color = cluster),
               size = 0.8, alpha = 0.7) +
    geom_line(data = yr_cpue,
              aes(x = doy, y = remaining_proportion_pct),
              color = "grey40", linewidth = 1.2, alpha = 0.6) +
    scale_color_manual(
      values = cluster_colors,
      breaks = cluster_order,
      labels = function(x) gsub("Cluster", "Group", x),
      guide  = guide_legend(override.aes = list(linewidth = 2))
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 0), "%"),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(MIN_DOY, MAX_DOY),
      expand = c(0, 0)
    ) +
    labs(
      title = paste("Cumulative distribution vs. Time —", yr),
      x     = "Date",
      y     = "Cumulative Percent of Group's Total Production",
      color = "Group"
    ) +
    base_theme

  fname <- file.path(PRES_DIR, paste0("cumulative_dist_with_cpue_", yr, ".png"))
  ggsave(fname, p, width = 7, height = 5, dpi = 300, bg = "white")
  cat("  Saved:", basename(fname), "\n")
}

cat("=== PRES-B COMPLETE:", length(available_years), "plots saved ===\n\n")

################################################################################
# DATA PREP FOR PRES-C  (replicates Analysis 4 & 6 logic exactly)
################################################################################

# ----- Analysis 4 cluster closure protection -----
closure_protection_by_unit <- cumulative_data %>%
  group_by(year, mgmt_river, cluster) %>%
  arrange(doy) %>%
  summarise(
    total_annual_production  = max(cumulative_production, na.rm = TRUE),
    closure_window_production = {
      closure_days <- doy[doy >= CLOSURE_START_DOY & doy <= CLOSURE_END_DOY]
      if (length(closure_days) > 0) {
        max_closure_doy <- max(closure_days)
        cumulative_production[doy == max_closure_doy][1]
      } else { 0 }
    },
    has_closure_data = any(doy >= CLOSURE_START_DOY & doy <= CLOSURE_END_DOY),
    .groups = "drop"
  ) %>%
  mutate(
    closure_window_production = ifelse(is.na(closure_window_production), 0, closure_window_production),
    unit_protection_pct = (closure_window_production / total_annual_production) * 100
  )

cluster_closure_protection <- closure_protection_by_unit %>%
  group_by(year, cluster) %>%
  summarise(
    cluster_total_production   = sum(total_annual_production,   na.rm = TRUE),
    cluster_closure_production = sum(closure_window_production, na.rm = TRUE),
    n_units = n(),
    .groups = "drop"
  ) %>%
  mutate(
    closure_protection_pct = (cluster_closure_production / cluster_total_production) * 100,
    cluster = factor(cluster, levels = rev(cluster_order))
  )

# ----- Analysis 5 offset protection -----
time_offsets <- list(
  "minus_3" = list(offset = -3, label = "-3 days (June 8)"),
  "plus_3"  = list(offset =  3, label = "+3 days (June 14)"),
  "plus_6"  = list(offset =  6, label = "+6 days (June 17)")
)

calculate_offset_protection <- function(cumulative_data, offset_days) {
  adj_start_doy <- CLOSURE_START_DOY + offset_days
  adj_end_doy   <- CLOSURE_END_DOY   + offset_days

  cumulative_data %>%
    group_by(year, mgmt_river, cluster) %>%
    arrange(doy) %>%
    summarise(
      total_annual_production   = max(cumulative_production, na.rm = TRUE),
      closure_window_production = {
        closure_days <- doy[doy >= adj_start_doy & doy <= adj_end_doy]
        if (length(closure_days) > 0) {
          max_closure_doy <- max(closure_days)
          cumulative_production[doy == max_closure_doy][1]
        } else { 0 }
      },
      has_closure_data = any(doy >= adj_start_doy & doy <= adj_end_doy),
      .groups = "drop"
    ) %>%
    mutate(
      closure_window_production = ifelse(is.na(closure_window_production), 0, closure_window_production),
      unit_protection_pct = (closure_window_production / total_annual_production) * 100
    ) %>%
    group_by(year, cluster) %>%
    summarise(
      cluster_total_production   = sum(total_annual_production,   na.rm = TRUE),
      cluster_closure_production = sum(closure_window_production, na.rm = TRUE),
      n_units = n(),
      .groups = "drop"
    ) %>%
    mutate(
      closure_protection_pct = (cluster_closure_production / cluster_total_production) * 100,
      cluster = factor(cluster, levels = rev(cluster_order))
    )
}

offset_results <- list()
for (offset_name in names(time_offsets)) {
  res <- calculate_offset_protection(cumulative_data, time_offsets[[offset_name]]$offset)
  offset_results[[offset_name]] <- res %>%
    mutate(offset_label = time_offsets[[offset_name]]$label)
}

combined_data <- bind_rows(
  cluster_closure_protection %>% mutate(offset_label = "Current (June 11)", offset_order = 2),
  offset_results$minus_3     %>% mutate(offset_order = 1),
  offset_results$plus_3      %>% mutate(offset_order = 3),
  offset_results$plus_6      %>% mutate(offset_order = 4)
) %>%
  mutate(offset_label = factor(offset_label,
                               levels = c("-3 days (June 8)", "Current (June 11)",
                                          "+3 days (June 14)", "+6 days (June 17)")))

# ----- CPUE foregone at each scenario end date -----
get_cpue_at_doy <- function(cpue_data, target_doy) {
  cpue_data %>%
    filter(doy == target_doy) %>%
    select(year, remaining_proportion_pct) %>%
    mutate(target_doy = target_doy,
           foregone_proportion_pct = 100 - remaining_proportion_pct)
}

cpue_scenarios <- bind_rows(
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY - 3) %>%
    mutate(offset_label = "-3 days (June 8)",   offset_order = 1, offset_short = "-3 days"),
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY)     %>%
    mutate(offset_label = "Current (June 11)",  offset_order = 2, offset_short = "0 days"),
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY + 3) %>%
    mutate(offset_label = "+3 days (June 14)",  offset_order = 3, offset_short = "+3 days"),
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY + 6) %>%
    mutate(offset_label = "+6 days (June 17)",  offset_order = 4, offset_short = "+6 days")
) %>%
  mutate(offset_label = factor(offset_label,
                               levels = c("-3 days (June 8)", "Current (June 11)",
                                          "+3 days (June 14)", "+6 days (June 17)")))

protection_cpue_combined <- combined_data %>%
  left_join(cpue_scenarios, by = c("year", "offset_label")) %>%
  filter(!is.na(foregone_proportion_pct)) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

mean_data <- protection_cpue_combined %>%
  group_by(cluster, offset_label) %>%
  summarise(
    n_years             = n(),
    mean_cpue_foregone  = mean(foregone_proportion_pct, na.rm = TRUE),
    sd_cpue_foregone    = sd(foregone_proportion_pct,   na.rm = TRUE),
    mean_protection     = mean(closure_protection_pct,  na.rm = TRUE),
    sd_protection       = sd(closure_protection_pct,    na.rm = TRUE),
    offset_short        = first(offset_short),
    .groups = "drop"
  ) %>%
  mutate(
    cluster      = factor(cluster, levels = cluster_order),
    offset_order = case_when(
      offset_label == "-3 days (June 8)"   ~ 1,
      offset_label == "Current (June 11)"  ~ 2,
      offset_label == "+3 days (June 14)"  ~ 3,
      offset_label == "+6 days (June 17)"  ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(cluster, offset_order)

################################################################################
# PRES-C: CPUE FOREGONE vs. PROTECTION — AXES FLIPPED
# x = Mean Protection Effectiveness (%)
# y = Mean Foregone Harvest Opportunity (%)
################################################################################

cat("=== PRES-C: CPUE foregone vs. protection (flipped axes) ===\n")

cluster_color_mapping <- setNames(cluster_colors, cluster_order)

p_flipped <- ggplot(mean_data,
                    aes(x    = mean_protection,
                        y    = mean_cpue_foregone,
                        color = cluster,
                        group = cluster)) +
  geom_line(linewidth = 1.5, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.6,
             stroke = 2, shape = 21,
             aes(fill = cluster), color = "black") +
  geom_text(data = mean_data %>%
              group_by(offset_label) %>%
              filter(mean_protection == max(mean_protection)),
            aes(label = offset_short),
            vjust = -1.2, hjust = 0.5, size = 6,
            show.legend = FALSE, fontface = "bold",
            color = "black") +
  scale_color_manual(values = cluster_color_mapping, name = "Group",
                     labels = paste0("Group ", 1:6)) +
  scale_fill_manual(values  = cluster_color_mapping, name = "Group",
                    labels  = paste0("Group ", 1:6)) +
  guides(color = guide_legend(title = "Group"),
         fill  = guide_legend(title = "Group")) +
  scale_x_continuous(
    name   = "Mean Protection Effectiveness (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = seq(0, 55, by = 10),
    limits = c(0, 55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = "Mean foregone harvest opportunity (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position      = "right",
    legend.title         = element_text(face = "bold", size = 20),
    legend.text          = element_text(face = "bold", size = 18),
    legend.background    = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin        = margin(6, 8, 6, 8),
    axis.title           = element_text(face = "bold", size = 20),
    axis.text            = element_text(face = "bold", size = 18),
    axis.line            = element_line(color = "black", linewidth = 1),
    panel.grid.major     = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor     = element_blank(),
    panel.border         = element_blank(),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(10, 20, 10, 10)
  )

fname_c <- file.path(PRES_DIR, "cpue_vs_protection_flipped_axes.png")
ggsave(fname_c, p_flipped, width = 12, height = 10, dpi = 300, bg = "white")
cat("  Saved:", basename(fname_c), "\n")
cat("=== PRES-C COMPLETE ===\n\n")

################################################################################
# SUMMARY
################################################################################

cat("================================\n")
cat("ALL PRESENTATION FIGURES COMPLETE\n")
cat("================================\n")
cat("Output directory:", PRES_DIR, "\n\n")
cat("PRES-A:", length(available_years), "individual cumulative plots (no CPUE)\n")
cat("  Files: cumulative_dist_no_cpue_<year>.png\n\n")
cat("PRES-B:", length(available_years), "individual cumulative plots (with CPUE)\n")
cat("  Files: cumulative_dist_with_cpue_<year>.png\n\n")
cat("PRES-C: 1 protection vs. foregone harvest plot (flipped axes)\n")
cat("  File: cpue_vs_protection_flipped_axes.png\n\n")
