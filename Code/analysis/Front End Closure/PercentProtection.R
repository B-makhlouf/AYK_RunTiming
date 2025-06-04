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

