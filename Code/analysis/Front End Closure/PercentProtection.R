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

