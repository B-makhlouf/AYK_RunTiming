# =============================================================================
# PROTECTION EFFECTIVENESS ANALYSIS
# Analyzes closure protection effectiveness and Q1 contribution patterns
# =============================================================================

# Load libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(sf)
library(tidyr)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

# Load closure summary data
closureSummary <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Frontend_Closure_Management/management_closure_summary_Kusko_2017_to_2021_06-01_to_06-11.csv")

# Load management river data
mgmt_data <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv") %>%
  filter(mgmt_river != "Johnson")  # Remove Johnson river

# Load spatial data
edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp", quiet = TRUE)
basin <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp", quiet = TRUE)

# =============================================================================
# PROTECTION EFFECTIVENESS CALCULATIONS
# =============================================================================

# Calculate protection effectiveness across harvest rate scenarios
harvest_rates <- c(0.1, 0.2, 0.3, 0.4, 0.5)  # 10% to 50%

protection_df <- closureSummary %>%
  crossing(harvest_rate = harvest_rates) %>%
  rename(
    year = Year,
    management_unit = Management_Unit,
    total_contribution_pct = Total_Contribution_Percent,
    percent_in_window = Percent_Of_Unit_In_Window
  ) %>%
  mutate(
    # Convert to proportions
    prop_in_window = percent_in_window / 100,
    prop_outside_window = 1 - prop_in_window,
    
    # Survival scenarios
    survival_no_closure = 1 - harvest_rate,
    survival_with_closure = prop_in_window + (prop_outside_window * (1 - harvest_rate)),
    
    # Protection metrics
    absolute_survival_benefit = survival_with_closure - survival_no_closure,
    relative_survival_improvement = (survival_with_closure - survival_no_closure) / survival_no_closure,
    percent_improvement = relative_survival_improvement * 100,
    protection_value = prop_in_window * harvest_rate
  )

# Filter for 30% harvest scenario (main analysis)
protection_30pct <- protection_df %>% filter(harvest_rate == 0.3)

# =============================================================================
# Q1 DATA PREPARATION
# =============================================================================

# Extract Q1 data
mgmt_q1_data <- mgmt_data %>% filter(quartile == "Q1")

# Calculate Q1 change (2018-2021) for units with ≥3% starting contribution
q1_change_data <- mgmt_q1_data %>%
  filter(year %in% c(2018, 2021)) %>%
  select(mgmt_river, year, within_quartile_prop) %>%
  pivot_wider(names_from = year, values_from = within_quartile_prop, names_prefix = "year_") %>%
  filter(year_2018 >= 0.02) %>%  # Filter to ≥2% starting contribution
  mutate(
    pct_change = ((year_2021 - year_2018) / year_2018) * 100,
    pct_change = ifelse(is.infinite(pct_change) | is.na(pct_change), 0, pct_change),
    starting_pct = year_2018 * 100,
    ending_pct = year_2021 * 100
  )


# =============================================================================
# SUMMARY DATA FOR SCATTER PLOTS
# =============================================================================

# Protection effectiveness summary
protection_summary <- protection_30pct %>%
  group_by(management_unit) %>%
  summarise(
    mean_protection = mean(percent_improvement, na.rm = TRUE),
    min_protection = min(percent_improvement, na.rm = TRUE),
    max_protection = max(percent_improvement, na.rm = TRUE),
    .groups = "drop"
  )

# Total run contribution summary
total_run_summary <- mgmt_data %>%
  group_by(mgmt_river) %>%
  summarise(mean_total_run = mean(total_run_prop * 100, na.rm = TRUE), .groups = "drop") %>%
  rename(management_unit = mgmt_river)

# Q1 contribution summary
q1_run_summary <- mgmt_q1_data %>%
  group_by(mgmt_river) %>%
  summarise(mean_q1_run = mean(within_quartile_prop * 100, na.rm = TRUE), .groups = "drop") %>%
  rename(management_unit = mgmt_river)

# Combined summaries
combined_summary <- protection_summary %>% left_join(total_run_summary, by = "management_unit")
combined_summary_q1 <- total_run_summary %>% left_join(q1_run_summary, by = "management_unit")

# =============================================================================
# COLOR PALETTES
# =============================================================================

# Create consistent color palettes
n_units <- length(unique(protection_30pct$management_unit))
unit_colors <- if (n_units <= 12) {
  brewer.pal(max(3, n_units), "Set3")
} else {
  colorRampPalette(brewer.pal(12, "Set3"))(n_units)
}

unit_colors_extended <- rainbow(n_units, s = 0.8, v = 0.8)
names(unit_colors_extended) <- sort(unique(protection_30pct$management_unit))

# Colors for stacked bar chart
stacked_colors <- if (length(unique(q1_stacked_data$mgmt_river)) <= 12) {
  brewer.pal(max(3, length(unique(q1_stacked_data$mgmt_river))), "Set3")
} else {
  colorRampPalette(brewer.pal(12, "Set3"))(length(unique(q1_stacked_data$mgmt_river)))
}
names(stacked_colors) <- levels(q1_stacked_data$mgmt_river)

# =============================================================================
# VISUALIZATION 1: PROTECTION EFFECTIVENESS TIME SERIES
# =============================================================================

p1 <- ggplot(protection_30pct, aes(x = year, y = percent_improvement, color = management_unit)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = unit_colors, name = "Management\nUnit") +
  scale_x_continuous(breaks = unique(protection_30pct$year)) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%"),
                     limits = c(0, max(protection_30pct$percent_improvement, na.rm = TRUE) * 1.05)) +
  labs(
    title = "Protection Effectiveness by Management Unit (30% Harvest Rate)",
    subtitle = "Percent improvement in survival due to June 1-11 closure window",
    x = "Year", y = "Percent Improvement in Survival",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "right", panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  )

print(p1)

# =============================================================================
# VISUALIZATION 2: Q1 CONTRIBUTION TIME SERIES
# =============================================================================

p2 <- ggplot(mgmt_q1_data, aes(x = year, y = within_quartile_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = unit_colors, name = "Management\nUnit") +
  scale_x_continuous(breaks = unique(mgmt_q1_data$year)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100, 1), "%"),
                     limits = c(0, max(mgmt_q1_data$within_quartile_prop, na.rm = TRUE) * 1.05)) +
  labs(
    title = "Q1 Proportional Contribution by Management Unit",
    subtitle = "Proportion of each management unit's production occurring in Q1",
    x = "Year", y = "Proportion of Unit's Production in Q1"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "right", panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  )

print(p2)

# =============================================================================
# VISUALIZATION 3: SPATIAL MAP OF Q1 CHANGE
# =============================================================================

# Join spatial data with change data
edges_with_change <- edges %>%
  left_join(q1_change_data, by = "mgmt_river") %>%
  filter(!is.na(mgmt_river) & mgmt_river != "")

p3 <- ggplot() +
  geom_sf(data = basin, fill = "gray95", color = "gray70", linewidth = 0.5, alpha = 0.3) +
  geom_sf(data = edges_with_change, aes(color = pct_change, linewidth = Str_Order), alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "% Change\n(2018-2021)",
    labels = function(x) paste0(round(x, 0), "%"),
    na.value = "gray60"
  ) +
  scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
  coord_sf(datum = NA) +
  labs(
    title = "Change in Q1 Contribution by Management Unit (2018-2021)",
    subtitle = "Blue = Decreased early-season contribution, Red = Increased early-season contribution",
    caption = "Based on proportional contribution within each management unit (≥2% starting contribution)"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front End Closure/Q1_contribution_change_map_2018_2021_filtered_3pct.png", 
       plot = p3, width = 12, height = 8, dpi = 300, bg = "white")



# =============================================================================
# VISUALIZATION 5: PROTECTION VS TOTAL RUN CONTRIBUTION
# =============================================================================

p4 <- ggplot(combined_summary, aes(x = mean_total_run, y = mean_protection, color = management_unit)) +
  geom_errorbar(aes(ymin = min_protection, ymax = max_protection), 
                width = 0.2, alpha = 0.7, linewidth = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = unit_colors_extended, name = "Management Unit") +
  scale_x_continuous(labels = function(x) paste0(round(x, 1), "%"),
                     limits = c(0, max(combined_summary$mean_total_run, na.rm = TRUE) * 1.1),
                     expand = expansion(mult = c(0.02, 0.1))) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  labs(
    title = "Mean Protection Effectiveness vs Mean Total Run Contribution",
    subtitle = "Error bars show min-max range of protection effectiveness across years (2017-2021)",
    x = "Mean Proportion of Total Run (%)", y = "Mean Protection Effectiveness (%)",
    caption = "Based on 30% harvest rate scenario and June 1-11 closure window"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom", legend.direction = "horizontal",
    panel.grid.minor = element_blank(), axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(color = guide_legend(ncol = 4))

ggsave("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front End Closure/mean_protection_vs_total_run_with_error_bars.png", 
       plot = p4, width = 12, height = 8, dpi = 300, bg = "white")

#####################

