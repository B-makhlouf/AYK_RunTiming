################################################################################
# 08_doy50_variability.R  --  INTERANNUAL VARIABILITY IN RUN TIMING BY GROUP
################################################################################
# STANDALONE ANALYSIS (additive; does not modify any other script).
#
# QUESTION
#   Are some regional stock groups more variable than others in WHEN their run
#   arrives? We quantify this with the standard deviation (SD, in days) of the
#   day-of-year (DOY) at which each group reaches 50% of its total seasonal run,
#   computed across all study years.
#
# METHOD (kept consistent with the rest of the pipeline)
#   1. Read the seasonal cumulative-distribution data written by
#      R/03_cumulative_distribution.R (config: CUMULATIVE_CSV).
#   2. Aggregate the management units into the six regional groups and build the
#      group-level cumulative curve using the SAME interval-based method as
#      Analysis 2 in R/06_cluster_analysis.R (so the curves are identical).
#   3. For each year x group, take the DOY-at-50% as the FIRST 3-day interval
#      whose group cumulative percent is >= 50% (matches the "days to 50%"
#      convention used in Analysis 3 of R/06).
#   4. Across years, summarise each group's DOY-at-50%: mean and SD (in days).
#      The group ranking by SD answers the question.
#
# NOTE ON THE METRIC
#   SD is reported in days. Because DOY is measured from Jan 1 (an arbitrary
#   origin), SD is origin-independent and is the natural, directly interpretable
#   measure of how much a group's 50%-run date shifts from year to year.
#
# OUTPUTS  (written to outputs/DOY50_Variability/ ; nothing else is touched)
#   CSV/doy50_by_year.csv             DOY-at-50% for every year x group
#   CSV/doy50_sd_summary.csv          per-group mean / SD, ranked most->least variable
#   doy50_variability_by_group.png    figure illustrating the result
################################################################################

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# ---- Configuration ----------------------------------------------------------
if (!exists("PROJECT_ROOT")) source("config.R")

CUMULATIVE_DATA_PATH <- CUMULATIVE_CSV        # written by R/03_cumulative_distribution.R

# Dedicated, self-contained output folder (kept separate from other analyses)
DOY50_DIR     <- file.path(OUTPUT_ROOT, "DOY50_Variability")
DOY50_CSV_DIR <- file.path(DOY50_DIR, "CSV")
dir.create(DOY50_CSV_DIR, recursive = TRUE, showWarnings = FALSE)

# Date window: identical to Analysis 2 in R/06 (June 1st - July 23rd)
MIN_DOY <- 149
MAX_DOY <- 205

# ---- Regional group (cluster) assignments -----------------------------------
# Identical mapping to R/06_cluster_analysis.R (copied here so this script is
# self-contained and does not depend on R/06 having been sourced).
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
cluster_order <- c("Cluster 1", "Cluster 2", "Cluster 3",
                   "Cluster 4", "Cluster 5", "Cluster 6")

cat("\n=== STARTING DOY-AT-50% VARIABILITY ANALYSIS ===\n")

# ---- Load cumulative distribution data --------------------------------------
cumulative_data <- read_csv(CUMULATIVE_DATA_PATH, show_col_types = FALSE) %>%
  filter(mgmt_river != "Johnson", mgmt_river != "Lower Kusko") %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster)) %>%
  filter(doy >= MIN_DOY & doy <= MAX_DOY)

cat("Loaded cumulative data; years:",
    paste(sort(unique(cumulative_data$year)), collapse = ", "), "\n")

# ---- Build group-level cumulative curve (same method as Analysis 2 in R/06) --
# Annual production total per group (sum of each unit's final cumulative value).
cluster_annual_totals <- cumulative_data %>%
  group_by(year, cluster, mgmt_river) %>%
  summarise(final_production = max(cumulative_production, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(year, cluster) %>%
  summarise(cluster_annual_total = sum(final_production, na.rm = TRUE),
            .groups = "drop")

# Interval-based aggregation -> group cumulative percent at each DOY.
cluster_aggregated <- cumulative_data %>%
  group_by(year, mgmt_river, cluster) %>%
  arrange(doy) %>%
  mutate(interval_production = cumulative_production -
           lag(cumulative_production, default = 0)) %>%
  ungroup() %>%
  group_by(year, doy, cluster) %>%
  summarise(cluster_interval_production = sum(interval_production, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(year, cluster) %>%
  arrange(doy) %>%
  mutate(cluster_cumulative_production = cumsum(cluster_interval_production)) %>%
  ungroup() %>%
  left_join(cluster_annual_totals, by = c("year", "cluster")) %>%
  mutate(cluster_cumulative_percent =
           (cluster_cumulative_production / cluster_annual_total) * 100)

# ---- DOY at 50% per year x group --------------------------------------------
# Convention: first 3-day interval whose group cumulative percent is >= 50%
# (mirrors the "days to 50%" rule used in Analysis 3 of R/06).
doy50_by_year <- cluster_aggregated %>%
  group_by(year, cluster) %>%
  arrange(doy) %>%
  summarise(
    doy_at_50 = {
      idx <- which(cluster_cumulative_percent >= 50)
      if (length(idx) > 0) doy[min(idx)] else NA_real_
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(doy_at_50)) %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  arrange(cluster, year)

# ---- SD across years, per group ---------------------------------------------
doy50_sd_summary <- doy50_by_year %>%
  group_by(cluster) %>%
  summarise(
    n_years    = n(),
    mean_doy50 = mean(doy_at_50),
    sd_doy50   = sd(doy_at_50),
    min_doy50  = min(doy_at_50),
    max_doy50  = max(doy_at_50),
    range_days = max(doy_at_50) - min(doy_at_50),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_doy50)) %>%                  # most -> least variable
  mutate(variability_rank = row_number())

# ---- Export CSVs ------------------------------------------------------------
write_csv(doy50_by_year,
          file.path(DOY50_CSV_DIR, "doy50_by_year.csv"))
write_csv(doy50_sd_summary,
          file.path(DOY50_CSV_DIR, "doy50_sd_summary.csv"))

# ---- Figure -----------------------------------------------------------------
# Horizontal range plot: each group's yearly DOY-at-50% (points), its across-year
# range (bar) and mean (white circle), with the SD labelled. Groups are ordered
# most -> least variable, so the figure reads top-to-bottom as the answer to
# "which groups are more variable in run timing?".
# Colour ramp matches the rest of the paper (R/06_cluster_analysis.R).
cluster_colors <- colorRampPalette(
  c("#8B0000", "#FF0000", "#FF8C00", "#87CEEB", "#1E90FF", "#000080"))(6)
names(cluster_colors) <- cluster_order

sd_order <- as.character(doy50_sd_summary$cluster)   # already sorted most -> least
group_label <- function(x) gsub("Cluster", "Group", x)

plot_pts <- doy50_by_year %>%
  mutate(cluster = factor(as.character(cluster), levels = sd_order))
plot_sum <- doy50_sd_summary %>%
  mutate(cluster  = factor(as.character(cluster), levels = sd_order),
         sd_label = sprintf("SD = %.1f days", sd_doy50))

x_right <- max(plot_pts$doy_at_50)

p_doy50 <- ggplot() +
  # across-year range
  geom_segment(data = plot_sum,
               aes(y = cluster, yend = cluster,
                   x = min_doy50, xend = max_doy50, color = cluster),
               linewidth = 5, alpha = 0.30, lineend = "round") +
  # individual years
  geom_point(data = plot_pts,
             aes(y = cluster, x = doy_at_50, color = cluster),
             size = 3, alpha = 0.85) +
  # across-year mean
  geom_point(data = plot_sum,
             aes(y = cluster, x = mean_doy50, fill = cluster),
             shape = 21, color = "white", size = 5, stroke = 1.6) +
  # SD label
  geom_text(data = plot_sum,
            aes(y = cluster, x = x_right + 1.5, label = sd_label),
            hjust = 0, size = 4, color = "gray25") +
  scale_color_manual(values = cluster_colors, guide = "none") +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_y_discrete(name = NULL, limits = rev(sd_order), labels = group_label) +
  scale_x_continuous(name = "Day of year at 50% of cumulative run",
                     breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.03, 0.20))) +
  labs(
    title = "Interannual variability in run timing by group",
    subtitle = paste0("Points = yearly DOY at 50% of run · white circle = mean ",
                      "· bar = across-year range\nGroups ordered most → least ",
                      "variable (standard deviation, SD, in days)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, color = "gray15"),
    plot.subtitle = element_text(size = 10.5, color = "gray35",
                                 margin = margin(b = 10)),
    axis.text.y   = element_text(size = 12, color = "gray20"),
    axis.text.x   = element_text(size = 11, color = "gray30"),
    axis.title.x  = element_text(size = 12, color = "gray20",
                                 face = "bold", margin = margin(t = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray92", linewidth = 0.4),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(16, 20, 16, 16)
  )

DOY50_FIG <- file.path(DOY50_DIR, "doy50_variability_by_group.png")
ggsave(DOY50_FIG, p_doy50, width = 9, height = 5.5, dpi = 300, bg = "white")

# ---- Console summary --------------------------------------------------------
cat("\n--- DOY-at-50% by group (SD across years) ---\n")
print(as.data.frame(
  doy50_sd_summary %>%
    mutate(across(c(mean_doy50, sd_doy50), ~round(.x, 3)))
), row.names = FALSE)

most  <- doy50_sd_summary$cluster[1]
least <- doy50_sd_summary$cluster[nrow(doy50_sd_summary)]
cat(sprintf("\nMost variable group:  %s (SD = %.2f days)\n",
            most,  doy50_sd_summary$sd_doy50[1]))
cat(sprintf("Least variable group: %s (SD = %.2f days)\n",
            least, doy50_sd_summary$sd_doy50[nrow(doy50_sd_summary)]))
cat("\nExported:\n  ", file.path(DOY50_CSV_DIR, "doy50_by_year.csv"),
    "\n  ", file.path(DOY50_CSV_DIR, "doy50_sd_summary.csv"),
    "\n  ", DOY50_FIG, "\n")
cat("=== DOY-AT-50% VARIABILITY ANALYSIS COMPLETE ===\n\n")
