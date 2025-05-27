################################################################################
# MAIN.R - SIMPLIFIED SALMON RUN TIMING ANALYSIS
################################################################################
# PURPOSE: Creates three types of maps:
#   1. DOY_Quartile: Shows proportion within each timing quartile
#   2. DOY_Total: Shows proportion of total annual run 
#   3. Front-end Closure: Shows protection during closure window
# NOTE: Removed all cumulative analysis and unnecessary complexity
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

# Source required files (keeping existing structure)
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))
source(here("code/doy_analysis.R"))
source(here("code/doy_total_analysis.R"))

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Clean previous output files
clean_outputs <- function() {
  dirs_to_clean <- c("DOY_Quartile", "DOY_Total", "FrontEnd_Closure")
  
  for (dir in dirs_to_clean) {
    full_path <- here("Basin Maps", dir)
    if (dir.exists(full_path)) {
      png_files <- list.files(full_path, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
      if (length(png_files) > 0) {
        file.remove(png_files)
        message(paste("Cleaned", length(png_files), "files from", dir))
      }
    }
  }
}

#' Get watershed-specific parameters
get_params <- function(watershed) {
  if (watershed == "Kusko") {
    list(sensitivity_threshold = 0.7, min_error = 0.0006, min_stream_order = 3)
  } else if (watershed == "Yukon") {
    list(sensitivity_threshold = 0.7, min_error = 0.003, min_stream_order = 5)
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon'")
  }
}

################################################################################
# ANALYSIS 1: DOY QUARTILE ANALYSIS
# PURPOSE: Shows what proportion of production occurs within each timing quartile
################################################################################

run_doy_analysis <- function(years, watersheds) {
  clean_outputs()
  
  for (watershed in watersheds) {
    for (year in years) {
      params <- get_params(watershed)
      
      message(paste("Processing DOY quartiles:", year, watershed))
      DOY_Quartile_Analysis(
        year = year,
        watershed = watershed,
        sensitivity_threshold = params$sensitivity_threshold,
        min_error = params$min_error,
        min_stream_order = params$min_stream_order,
        HUC = 8
      )
    }
  }
}

################################################################################
# ANALYSIS 2: DOY TOTAL ANALYSIS  
# PURPOSE: Shows what proportion of the TOTAL ANNUAL RUN each quartile represents
################################################################################

run_doy_total_analysis <- function(years, watersheds) {
  for (watershed in watersheds) {
    for (year in years) {
      params <- get_params(watershed)
      
      message(paste("Processing DOY total:", year, watershed))
      DOY_Total_Analysis(
        year = year,
        watershed = watershed,
        sensitivity_threshold = params$sensitivity_threshold,
        min_error = params$min_error,
        min_stream_order = params$min_stream_order,
        HUC = 8
      )
    }
  }
}

################################################################################
# ANALYSIS 3: FRONT-END CLOSURE ANALYSIS
# PURPOSE: Shows what proportion of run is protected during closure window
################################################################################

run_closure_analysis <- function(years, watershed = "Kusko", 
                                 closure_start = "06-01", closure_end = "06-11") {
  
  # Use existing closure function from FrontEndClosure.R
  source(here("code/FrontEndClosure.R"))
  
  analyze_closure_window(
    years = years,
    watershed = watershed,
    closure_start = closure_start,
    closure_end = closure_end,
    create_maps = TRUE
  )
}

################################################################################
# MAIN EXECUTION SECTION
################################################################################

# Set analysis parameters
years <- c("2017","2018", "2019", "2020", "2021")
watersheds <- c("Kusko")

message("=== Running Simplified Analysis Suite ===")

# 1. DOY Quartiles (proportion within each quartile)
message("Starting DOY Quartile Analysis...")
run_doy_analysis(years, watersheds)

# 2. DOY Total (proportion of total run) 
message("Starting DOY Total Analysis...")
run_doy_total_analysis(years, watersheds)

# 3. Front-end Closure
message("Starting Front-end Closure Analysis...")
run_closure_analysis(years, "Kusko")

message("=== Analysis Complete ===")