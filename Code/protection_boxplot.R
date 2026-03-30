################################################################################
# FRONT-END CLOSURE PROTECTION — BOXPLOT BY GROUP
# Uses same data pipeline as PresentationFigures.R (PRES-D section)
# X-axis: Group (Group 1–6)
# Y-axis: Front-End Closure Protection (%)
# Data filtered to Current (June 11) offset only
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

################################################################################
# SETUP — mirrors PresentationFigures.R exactly
################################################################################

OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming"
PRES_DIR   <- file.path(OUTPUT_DIR, "Figures", "PresFigures")
CSV_DIR    <- file.path(OUTPUT_DIR, "CSV_Exports")

CUMULATIVE_DATA_PATH <- file.path(OUTPUT_DIR,
  "Analysis_Results/Cumulative_Distribution/CSV/cumulative_distribution_data_Kusko_3day_intervals.csv")

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

cluster_order  <- c("Cluster 1","Cluster 2","Cluster 3",
                    "Cluster 4","Cluster 5","Cluster 6")
cluster_labels <- paste0("Group ", 1:6)
names(cluster_labels) <- cluster_order

cluster_colors <- colorRampPalette(
  c("#8B0000","#FF0000","#FF8C00","#87CEEB","#1E90FF","#000080")
)(6)
names(cluster_colors) <- cluster_order

################################################################################
# DATA LOADING
################################################################################

MIN_DOY <- 149
MAX_DOY <- 205

cumulative_data <- read_csv(CUMULATIVE_DATA_PATH) %>%
  filter(mgmt_river != "Johnson", mgmt_river != "Lower Kusko") %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster)) %>%
  filter(doy >= MIN_DOY & doy <= MAX_DOY)

################################################################################
# CALCULATE CLOSURE PROTECTION (current offset only)
################################################################################

closure_protection_by_unit <- cumulative_data %>%
  group_by(year, mgmt_river, cluster) %>%
  arrange(doy) %>%
  summarise(
    total_annual_production   = max(cumulative_production, na.rm = TRUE),
    closure_window_production = {
      closure_days <- doy[doy >= CLOSURE_START_DOY & doy <= CLOSURE_END_DOY]
      if (length(closure_days) > 0) {
        max_closure_doy <- max(closure_days)
        cumulative_production[doy == max_closure_doy][1]
      } else { 0 }
    },
    .groups = "drop"
  ) %>%
  mutate(
    closure_window_production = ifelse(is.na(closure_window_production), 0,
                                       closure_window_production)
  )

cluster_closure_protection <- closure_protection_by_unit %>%
  group_by(year, cluster) %>%
  summarise(
    cluster_total_production   = sum(total_annual_production,   na.rm = TRUE),
    cluster_closure_production = sum(closure_window_production, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    closure_protection_pct = (cluster_closure_production / cluster_total_production) * 100,
    cluster = factor(cluster, levels = cluster_order)
  )

################################################################################
# BOXPLOT
################################################################################

p_box <- ggplot(cluster_closure_protection,
                aes(x = cluster, y = closure_protection_pct,
                    fill = cluster, color = cluster)) +
  geom_boxplot(
    width        = 0.55,
    alpha        = 0.55,
    outlier.shape = NA,   # hide default outlier points; jitter shows all
    linewidth    = 0.8
  ) +
  geom_jitter(
    width  = 0.15,
    size   = 2.5,
    alpha  = 0.6,
    shape  = 21,
    stroke = 0.6,
    aes(fill = cluster),
    color  = "white"
  ) +
  scale_fill_manual(values  = cluster_colors, guide = "none") +
  scale_color_manual(values = cluster_colors, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  scale_y_continuous(
    name   = "Front-End Closure Protection (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Front-End Closure Protection by Group",
    x     = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title         = element_text(face = "bold", size = 13, color = "gray20"),
    axis.text          = element_text(size = 12, color = "gray20"),
    axis.line.x        = element_line(color = "gray40", linewidth = 0.5),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(10, 20, 10, 10)
  )

################################################################################
# SAVE
################################################################################

out_path <- file.path(PRES_DIR, "protection_boxplot_by_group.png")
ggsave(out_path, p_box, width = 11, height = 7, dpi = 300, bg = "white")
cat("Saved:", out_path, "\n")
