################################################################################
# MAIN.R - SALMON RUN TIMING ANALYSIS SUITE (CLEANED - NO DFA)
################################################################################
# PURPOSE: Creates management river maps and exports CSV data:
#   1. DOY_Quartile: Shows proportion within each timing quartile
#   2. Front-end Closure: Shows protection during closure window
#   3. Average Quartile Maps: Shows average production across all years
#   4. Cumulative Distribution: Shows timing progression by management unit
# NOTE: DFA-related functions removed - only core spatial analysis
################################################################################

################################################################################
# SETUP: LOAD LIBRARIES AND SOURCE FILES
################################################################################
library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyr)
library(scales)
library(readr)
library(png)

# Source required files for spatial analysis only
source(here("Code/utils/spatial_utils.R"))
source(here("Code/utils/visualization.R"))
source(here("Code/assignment.R"))
source(here("Code/doy_analysis.R"))

################################################################################
# CORRECTED WATERSHED ORDERING FUNCTIONS
################################################################################

#' Get standardized watershed ordering for all plots and analyses (CORRECTED)
get_watershed_order <- function() {
  # Exact order from your document - upstream to downstream
  watershed_order <- c(
    "N. Fork Kusko",
    "E. Fork Kuskokwim River", 
    "S. Fork Kusko",
    "Upper Kusko Main",
    "Big River",
    "Takotna and Nixon Fork",
    "Tatlawiksuk",
    "Swift",
    "Stony", 
    "Holitna",
    "Hoholitna",
    "Middle Kusko Main",
    "George",
    "Oskakawlik", 
    "Holokuk",
    "Aniak",
    "Tuluksak",
    "Kisaralik",
    "Kwethluk",
    "Johnson",
    "Lower Kusko"
  )
  
  return(watershed_order)
}

################################################################################
# MAIN EXECUTION SECTION
################################################################################

# Set analysis parameters
years <- c(2017, 2018, 2019, 2020, 2021, 2022)
watersheds <- c("Kusko")

message("=== Running Core Salmon Run Timing Analysis Suite (CORRECTED WATERSHED ORDERING) ===")

# 1. DOY Quartiles (Management Rivers only) + CSV export
run_doy_analysis(years, watersheds, export_csv = TRUE)

# 2. Front-end Closure Boxplot Analysis (CORRECTED ORDERING)
message("2. Starting Front-end Closure Boxplot Analysis (CORRECTED ORDERING)...")
run_closure_boxplot_analysis(years, "Kusko")

# 3. Average Quartile Maps (CORRECTED ORDERING)
message("3. Starting Average Quartile Maps Analysis (CORRECTED ORDERING)...")
run_average_quartile_analysis(years, "Kusko")

# 4. Cumulative Distribution Analysis
message("4. Starting Cumulative Distribution Analysis...")
run_cumulative_distribution_analysis(years, "Kusko", interval_days = 3)

