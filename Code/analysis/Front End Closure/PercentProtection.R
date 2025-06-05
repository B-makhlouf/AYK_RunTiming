# Protection Effectiveness Calculation - Simplified Data Frame
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Load your closure summary data
closureSummary <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Frontend_Closure_Management/management_closure_summary_Kusko_2017_to_2021_06-01_to_06-11.csv")

# Define harvest rates to analyze
harvest_rates <- c(0.1, 0.2, 0.3, 0.4, 0.5)  # 10% to 50%

# Create the protection effectiveness data frame
protection_df <- closureSummary %>%
  # Expand to include all harvest rate scenarios
  crossing(harvest_rate = harvest_rates) %>%
  # Clean column names
  rename(
    year = Year,
    management_unit = Management_Unit,
    total_contribution_pct = Total_Contribution_Percent,
    percent_in_window = Percent_Of_Unit_In_Window
  ) %>%
  mutate(
    # Create scenario label
    harvest_scenario = paste0(harvest_rate * 100, "%"),
    
    # Convert to proportions
    prop_in_window = percent_in_window / 100,
    prop_outside_window = 1 - prop_in_window,
    
    # Scenario 1: No closure - harvest applies to all fish
    survival_no_closure = 1 - harvest_rate,
    
    # Scenario 2: With closure - fish in window protected, rest face harvest
    survival_with_closure = prop_in_window + (prop_outside_window * (1 - harvest_rate)),
    
    # Protection effectiveness metrics
    absolute_survival_benefit = survival_with_closure - survival_no_closure,
    relative_survival_improvement = (survival_with_closure - survival_no_closure) / survival_no_closure,
    
    # Additional useful metrics
    percent_improvement = relative_survival_improvement * 100,
    protection_value = prop_in_window * harvest_rate  # Direct protection value
  ) %>%
  # Reorder columns for clarity
  select(
    year, management_unit, harvest_scenario, harvest_rate,
    total_contribution_pct, percent_in_window, prop_in_window,
    survival_no_closure, survival_with_closure,
    absolute_survival_benefit, relative_survival_improvement, 
    percent_improvement, protection_value
  ) %>%
  arrange(year, management_unit, harvest_rate)

###################################################
######### Figure 1. Ts protection for 30% by management unit 

# Filter for 30% harvest scenario
protection_30pct <- protection_df %>%
  filter(harvest_rate == 0.3)  # 30% harvest rate

# Create color palette for management units
n_units <- length(unique(protection_30pct$management_unit))
unit_colors <- if (n_units <= 12) {
  brewer.pal(max(3, n_units), "Set3")
} else {
  colorRampPalette(brewer.pal(12, "Set3"))(n_units)
}

# Create the time series plot
p1 <- ggplot(protection_30pct, aes(x = year, y = percent_improvement, color = management_unit)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = unit_colors, name = "Management\nUnit") +
  scale_x_continuous(breaks = unique(protection_30pct$year)) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, max(protection_30pct$percent_improvement, na.rm = TRUE) * 1.05)
  ) +
  labs(
    title = "Protection Effectiveness by Management Unit (30% Harvest Rate)",
    subtitle = "Percent improvement in survival due to June 1-11 closure window",
    x = "Year",
    y = "Percent Improvement in Survival",
    caption = "Higher values = closure more effective for that management unit"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  )

print(p1)

###################################################
######### Figure 2. Q1 Management Unit Contribution Over Time
# Filter for Q1 data only
mgmt_q1_data <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv") %>%
  filter(quartile == "Q1") %>%
  filter(mgmt_river != "Johnson")  # Remove Johnson as in your protection code

# Create Q1 time series plot
p2 <- ggplot(mgmt_q1_data, aes(x = year, y = within_quartile_prop, color = mgmt_river)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = unit_colors, name = "Management\nUnit") +
  scale_x_continuous(breaks = unique(mgmt_q1_data$year)) +
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100, 1), "%"),
    limits = c(0, max(mgmt_q1_data$within_quartile_prop, na.rm = TRUE) * 1.05)
  ) +
  labs(
    title = "Q1 Proportional Contribution by Management Unit",
    subtitle = "Proportion of each management unit's production occurring in Q1",
    x = "Year",
    y = "Proportion of Unit's Production in Q1",
    caption = "Higher values = management unit contributes more to early-season run timing"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  )

print(p2)

#################################################################################
###################################################
######### Figure 3. Map of % Change in Q1 Contribution (2017-2021)
###################################################
######### Figure 3. Map of % Change in Q1 Contribution (2017-2021)

# Calculate % change in Q1 contribution from 2017 to 2021
q1_change_data <- mgmt_q1_data %>%
  filter(year %in% c(2017, 2021)) %>%
  select(mgmt_river, year, within_quartile_prop) %>%
  pivot_wider(names_from = year, values_from = within_quartile_prop, names_prefix = "year_") %>%
  mutate(
    pct_change = ((year_2021 - year_2017) / year_2017) * 100,
    pct_change = ifelse(is.infinite(pct_change) | is.na(pct_change), 0, pct_change)
  ) %>%
  select(mgmt_river, pct_change)

# Load spatial data
edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp", quiet = TRUE)
basin <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp", quiet = TRUE)

# Join change data with spatial data
edges_with_change <- edges %>%
  left_join(q1_change_data, by = "mgmt_river") %>%
  filter(!is.na(mgmt_river) & mgmt_river != "")

# Create the map
p3 <- ggplot() +
  geom_sf(data = basin, fill = "gray95", color = "gray70", linewidth = 0.5, alpha = 0.3) +
  geom_sf(data = edges_with_change, aes(color = pct_change, linewidth = Str_Order), alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "% Change\n(2017-2021)",
    labels = function(x) paste0(round(x, 0), "%"),
    na.value = "gray60"
  ) +
  scale_linewidth_continuous(
    range = c(0.5, 2.5), 
    name = "Stream\nOrder",
    guide = "none"  # Hide stream order legend to keep it simple
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Change in Q1 Contribution by Management Unit (2017-2021)",
    subtitle = "Blue = Decreased early-season contribution, Red = Increased early-season contribution",
    caption = "Based on proportional contribution within each management unit"
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

# Save as PNG
ggsave("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front End Closure/Q1_contribution_change_map_2017_2021.png", 
       plot = p3, width = 12, height = 8, dpi = 300, bg = "white")

#######################
###################################################
######### Figure 4. Mean Protection vs Mean Total Run Proportion (with error bars)

# Calculate mean protection effectiveness (keep this the same)
protection_summary <- protection_30pct %>%
  group_by(management_unit) %>%
  summarise(
    mean_protection = mean(percent_improvement, na.rm = TRUE),
    min_protection = min(percent_improvement, na.rm = TRUE),
    max_protection = max(percent_improvement, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate mean total run contribution from the management river data
total_run_summary <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv") %>%
  filter(mgmt_river != "Johnson") %>%
  group_by(mgmt_river) %>%
  summarise(
    mean_total_run = mean(total_run_prop * 100, na.rm = TRUE),  # Convert to percentage
    .groups = "drop"
  ) %>%
  rename(management_unit = mgmt_river)

# Combine the summaries
combined_summary <- protection_summary %>%
  left_join(total_run_summary, by = "management_unit")

# Create a better color palette for many classes
n_units <- length(unique(combined_summary$management_unit))
unit_colors_extended <- rainbow(n_units, s = 0.8, v = 0.8)
names(unit_colors_extended) <- sort(unique(combined_summary$management_unit))

# Create scatter plot with error bars only on protection effectiveness
p4 <- ggplot(combined_summary, aes(x = mean_total_run, y = mean_protection, color = management_unit)) +
  geom_errorbar(aes(ymin = min_protection, ymax = max_protection), 
                width = 0.2, alpha = 0.7, linewidth = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = unit_colors_extended, name = "Management Unit") +
  scale_x_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, max(combined_summary$mean_total_run, na.rm = TRUE) * 1.1),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  labs(
    title = "Mean Protection Effectiveness vs Mean Total Run Contribution",
    subtitle = "Error bars show min-max range of protection effectiveness across years (2017-2021)",
    x = "Mean Proportion of Total Run (%)",
    y = "Mean Protection Effectiveness (%)",
    caption = "Based on 30% harvest rate scenario and June 1-11 closure window"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(color = guide_legend(ncol = 4))  # Arrange legend in 4 columns

# Save as PNG
ggsave("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front End Closure/mean_protection_vs_total_run_with_error_bars.png", 
       plot = p4, width = 12, height = 8, dpi = 300, bg = "white")

###################################################
######### Figure 5. Total Run Contribution vs Q1 Contribution

# Calculate mean total run contribution (sum across all quartiles for each year, then average)
total_run_summary <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv") %>%
  filter(mgmt_river != "Johnson") %>%
  group_by(mgmt_river, year) %>%
  summarise(
    yearly_total_run = sum(total_run_prop, na.rm = TRUE) * 100,  # Sum across quartiles, convert to %
    .groups = "drop"
  ) %>%
  group_by(mgmt_river) %>%
  summarise(
    mean_total_run = mean(yearly_total_run, na.rm = TRUE),  # Average across years
    .groups = "drop"
  ) %>%
  rename(management_unit = mgmt_river)

# Calculate mean Q1 contribution using within_quartile_prop (proportion of Q1 run)
q1_run_summary <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv") %>%
  filter(mgmt_river != "Johnson") %>%
  filter(quartile == "Q1") %>%
  group_by(mgmt_river) %>%
  summarise(
    mean_q1_run = mean(within_quartile_prop * 100, na.rm = TRUE),  # This should sum to 100%
    .groups = "drop"
  ) %>%
  rename(management_unit = mgmt_river)

# Combine the summaries
combined_summary_q1 <- total_run_summary %>%
  left_join(q1_run_summary, by = "management_unit")

# Create a better color palette for many classes
n_units <- length(unique(combined_summary_q1$management_unit))
unit_colors_extended <- rainbow(n_units, s = 0.8, v = 0.8)
names(unit_colors_extended) <- sort(unique(combined_summary_q1$management_unit))

# Create scatter plot 
p5 <- ggplot(combined_summary_q1, aes(x = mean_total_run, y = mean_q1_run, color = management_unit)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_color_manual(values = unit_colors_extended, name = "Management Unit") +
  scale_x_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, max(combined_summary_q1$mean_total_run, na.rm = TRUE) * 1.1),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  labs(
    title = "Total Run Contribution vs Q1 Contribution by Management Unit",
    subtitle = "Dashed line shows equal contribution (Q1 = Total)",
    x = "Mean Proportion of Total Run (%)",
    y = "Mean Proportion of Q1 Run (%)",
    caption = "Points above line = higher Q1 contribution, points below = lower Q1 contribution"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(color = guide_legend(ncol = 4))  # Arrange legend in 4 columns

# Save as PNG
ggsave("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front End Closure/total_run_vs_q1_contribution.png", 
       plot = p5, width = 12, height = 8, dpi = 300, bg = "white")