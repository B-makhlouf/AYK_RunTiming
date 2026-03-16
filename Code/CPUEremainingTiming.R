################################################################################
# SIMPLIFIED REMAINING CPUE CALCULATION
################################################################################
# PURPOSE: Plot daily CPUE and calculate remaining proportion at each day
# OUTPUT: Single CSV with all data, faceted CPUE curve figure with DOY
################################################################################

library(dplyr)
library(readr)
library(ggplot2)

################################################################################
# PARAMETERS
################################################################################

years <- c(2017, 2018, 2019, 2020, 2021, 2022)
watershed <- "Kusko"

# Output directory
output_dir <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming/Analysis_Results/RemainingCPUE"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# LOAD DATA
################################################################################

cat("\n=== LOADING CPUE DATA ===\n\n")

all_years_data <- list()

for (year in years) {
  file_path <- paste0(
    "/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/",
    "Arctic_Yukon_Kuskokwim_Data/Data/Natal Origin Analysis Data/",
    "03_Natal Origins Genetics CPUE/",
    year, "_", watershed, "_Natal_Origins_Genetics_CPUE.csv"
  )
  
  if (file.exists(file_path)) {
    # Load the data
    cpue <- read.csv(file_path)
    
    # Get unique combinations of DOY and dailyCPUEprop (aggregate by day)
    # This handles the fact that there are multiple fish per day
    cpue_daily <- cpue %>%
      select(DOY, dailyCPUEprop) %>%
      distinct() %>%
      rename(doy = DOY, cpue = dailyCPUEprop) %>%
      arrange(doy) %>%
      mutate(year = year)
    
    # Replace any NA with a 0 
    cpue_daily$cpue[is.na(cpue_daily$cpue)] <- 0
    
    all_years_data[[as.character(year)]] <- cpue_daily
    cat("   ✓ Loaded year", year, "-", nrow(cpue_daily), "unique days\n")
  } else {
    cat("   ✗ File not found for year", year, "\n")
  }
}

cpue_data <- bind_rows(all_years_data)
cat("\n   Total daily records:", nrow(cpue_data), "\n\n")

################################################################################
# CALCULATE REMAINING PROPORTION
################################################################################

cat("=== CALCULATING REMAINING PROPORTION ===\n\n")

cpue_analysis <- cpue_data %>%
  arrange(year, doy) %>%
  group_by(year) %>%
  mutate(
    total_annual_cpue = sum(cpue, na.rm = TRUE),
    cumulative_cpue = cumsum(cpue),
    cumulative_proportion = cumulative_cpue / total_annual_cpue,
    remaining_cpue = total_annual_cpue - cumulative_cpue,
    remaining_proportion = 1 - cumulative_proportion
  ) %>%
  ungroup() %>%
  select(year, doy, cpue, total_annual_cpue, cumulative_cpue, 
         cumulative_proportion, remaining_cpue, remaining_proportion)

cat("   ✓ Calculated remaining proportion for each day\n\n")

################################################################################
# EXPORT SINGLE CSV
################################################################################

output_file <- file.path(output_dir, "daily_cpue_with_remaining_proportion.csv")
write_csv(cpue_analysis, output_file)

cat("=== SAVED OUTPUT ===\n")
cat("   ✓ File:", output_file, "\n\n")

################################################################################
# FACETED CPUE CURVE FIGURE - SINGLE COLUMN WITH DOY
################################################################################

cat("=== GENERATING FIGURE ===\n\n")

p <- ggplot(cpue_analysis, aes(x = doy, y = cpue)) +
  geom_line(size = 1, color = "steelblue") +
  facet_wrap(~ year, ncol = 1, scales = "fixed") +
  labs(
    title = "Daily CPUE by Year",
    x = "Day of Year",
    y = "Daily CPUE Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(output_dir, "cpue_curves_by_year_faceted.png"),
       p, width = 10, height = 12, dpi = 300, bg = "white")

cat("   ✓ Saved cpue_curves_by_year_faceted.png\n\n")

################################################################################
# SUMMARY
################################################################################

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("=== COMPLETE ===\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
cat("Years analyzed:", paste(unique(cpue_analysis$year), collapse = ", "), "\n")
cat("Days per year:\n")
print(cpue_analysis %>% count(year))
cat("\nOutputs saved to:", output_dir, "\n\n")