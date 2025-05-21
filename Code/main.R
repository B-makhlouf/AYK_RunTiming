# main.R
# Main execution script for DOY timing analysis

library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(tidyverse)

# Source all required function files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))
source(here("code/doy_analysis.R"))
source(here("code/cumulative_doy_analysis.R"))
source(here("code/export_doy_values.R"))
source(here("code/doy_total_analysis.R"))
source(here("code/cumulative_doy_total_analysis.R"))
source(here("code/export_doy_total_values.R"))

#' Clean previous DOY output files
#'
#' @param base_dir Base directory to clean
#' @return Number of files removed
clean_doy_outputs <- function(base_dir = here("Basin Maps")) {
  message("Cleaning previous DOY outputs...")
  
  # Define directories to clean
  doy_dirs <- c(
    "DOY_Quartile",
    "DOY_Cumulative",
    "DOY_Total",
    "DOY_Cumulative_Total"
  )
  
  total_removed <- 0
  for (dir in doy_dirs) {
    dir_path <- file.path(base_dir, dir)
    if (dir.exists(dir_path)) {
      png_files <- list.files(dir_path, 
                              pattern = "\\.png$", 
                              recursive = TRUE,
                              full.names = TRUE)
      
      if (length(png_files) > 0) {
        file.remove(png_files)
        total_removed <- total_removed + length(png_files)
      }
    }
  }
  
  message(paste(total_removed, "PNG files removed."))
  return(total_removed)
}

#' Get watershed-specific parameters
#'
#' @param watershed Character: "Kusko" or "Yukon"
#' @return List of parameters for the watershed
get_watershed_params <- function(watershed) {
  if (watershed == "Yukon") {
    return(list(
      sensitivity_threshold = 0.7,
      min_error = 0.003,
      min_stream_order = 5
    ))
  } else if (watershed == "Kusko") {
    return(list(
      sensitivity_threshold = 0.7,
      min_error = 0.0006,
      min_stream_order = 3
    ))
  } else {
    stop(paste("Unknown watershed:", watershed))
  }
}

#' Process a single dataset for DOY timing analysis
#'
#' @param dataset Dataset name in the format "year_watershed"
#' @param run_cumulative Whether to run cumulative DOY analysis
#' @return Invisibly returns NULL
process_doy_dataset <- function(dataset, run_cumulative = TRUE) {
  # Extract year and watershed
  parts <- strsplit(dataset, "_")[[1]]
  year <- parts[1]
  watershed <- parts[2]
  
  message(paste("Processing DOY analysis for:", year, watershed))
  
  # Get watershed-specific parameters
  params <- get_watershed_params(watershed)
  
  # Run DOY quartile analysis
  message("  Processing DOY quartiles")
  DOY_Quartile_Analysis(
    year = year,
    watershed = watershed,
    sensitivity_threshold = params$sensitivity_threshold,
    min_error = params$min_error,
    min_stream_order = params$min_stream_order,
    HUC = 8
  )
  
  # Run cumulative DOY analysis if requested
  if (run_cumulative) {
    message("  Processing cumulative DOY analysis")
    Cumulative_DOY_Analysis(
      year = year,
      watershed = watershed,
      sensitivity_threshold = params$sensitivity_threshold,
      min_error = params$min_error,
      min_stream_order = params$min_stream_order,
      HUC = 8
    )
  }
  
  # Clean up any variables that might cause conflicts in subsequent runs
  gc() # Force garbage collection
  
  return(invisible(NULL))
}

#' Process a single dataset for DOY total timing analysis
#'
#' @param dataset Dataset name in the format "year_watershed"
#' @param run_cumulative Whether to run cumulative DOY analysis
#' @return Invisibly returns NULL
process_doy_total_dataset <- function(dataset, run_cumulative = TRUE) {
  # Extract year and watershed
  parts <- strsplit(dataset, "_")[[1]]
  year <- parts[1]
  watershed <- parts[2]
  
  message(paste("Processing DOY total analysis for:", year, watershed))
  
  # Get watershed-specific parameters
  params <- get_watershed_params(watershed)
  
  # Run DOY total quartile analysis
  message("  Processing DOY total quartiles")
  DOY_Total_Analysis(
    year = year,
    watershed = watershed,
    sensitivity_threshold = params$sensitivity_threshold,
    min_error = params$min_error,
    min_stream_order = params$min_stream_order,
    HUC = 8
  )
  
  # Run cumulative DOY total analysis if requested
  if (run_cumulative) {
    message("  Processing cumulative DOY total analysis")
    Cumulative_DOY_Total_Analysis(
      year = year,
      watershed = watershed,
      sensitivity_threshold = params$sensitivity_threshold,
      min_error = params$min_error,
      min_stream_order = params$min_stream_order,
      HUC = 8
    )
  }
  
  # Clean up any variables that might cause conflicts in subsequent runs
  gc() # Force garbage collection
  
  return(invisible(NULL))
}

#' Run DOY timing analysis on specific datasets
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param run_cumulative Whether to run cumulative DOY analysis
#' @return Invisibly returns NULL
run_doy_timing_analysis <- function(years, watersheds, run_cumulative = TRUE) {
  # Create all combinations of years and watersheds
  datasets <- c()
  for (year in years) {
    for (watershed in watersheds) {
      datasets <- c(datasets, paste(year, watershed, sep = "_"))
    }
  }
  
  message(paste("Processing DOY timing analysis for datasets:"))
  message(paste(datasets, collapse = ", "))
  
  # Clean outputs
  clean_doy_outputs()
  
  # Process each dataset
  for (dataset in datasets) {
    tryCatch({
      process_doy_dataset(dataset, run_cumulative)
      message(paste("Successfully processed DOY analysis for:", dataset))
    }, error = function(e) {
      message(paste("Error processing DOY analysis for", dataset, ":", e$message))
    })
  }
  
  message("DOY timing analysis complete!")
  return(invisible(NULL))
}

#' Run DOY total timing analysis on specific datasets
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param run_cumulative Whether to run cumulative DOY analysis
#' @return Invisibly returns NULL
run_doy_total_timing_analysis <- function(years, watersheds, run_cumulative = TRUE) {
  # Create all combinations of years and watersheds
  datasets <- c()
  for (year in years) {
    for (watershed in watersheds) {
      datasets <- c(datasets, paste(year, watershed, sep = "_"))
    }
  }
  
  message(paste("Processing DOY total timing analysis for datasets:"))
  message(paste(datasets, collapse = ", "))
  
  # Clean outputs (optional, already done in run_doy_timing_analysis if run together)
  # clean_doy_outputs()
  
  # Process each dataset
  for (dataset in datasets) {
    tryCatch({
      process_doy_total_dataset(dataset, run_cumulative)
      message(paste("Successfully processed DOY total analysis for:", dataset))
    }, error = function(e) {
      message(paste("Error processing DOY total analysis for", dataset, ":", e$message))
    })
  }
  
  message("DOY total timing analysis complete!")
  return(invisible(NULL))
}

#' Run DOY timing analysis and export results
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param run_cumulative Whether to run cumulative analysis
#' @param export_data Whether to export data to CSV
#' @return Invisibly returns NULL
run_full_doy_analysis <- function(years, watersheds, 
                                  run_cumulative = TRUE,
                                  export_data = TRUE) {
  # Run the timing analysis
  run_doy_timing_analysis(years = years, 
                          watersheds = watersheds,
                          run_cumulative = run_cumulative)
  
  # Export the results if requested
  if (export_data) {
    message("Exporting DOY analysis results...")
    export_doy_analysis_results(
      years = years,
      watersheds = watersheds,
      include_cumulative = run_cumulative
    )
  }
  
  message("Full DOY analysis and export complete!")
  return(invisible(NULL))
}

#' Run DOY total timing analysis and export results
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param run_cumulative Whether to run cumulative analysis
#' @param export_data Whether to export data to CSV
#' @return Invisibly returns NULL
run_full_doy_total_analysis <- function(years, watersheds, 
                                        run_cumulative = TRUE,
                                        export_data = TRUE) {
  # Run the timing analysis
  run_doy_total_timing_analysis(years = years, 
                                watersheds = watersheds,
                                run_cumulative = run_cumulative)
  
  # Export the results if requested
  if (export_data) {
    message("Exporting DOY total analysis results...")
    export_doy_total_analysis_results(
      years = years,
      watersheds = watersheds,
      include_cumulative = run_cumulative
    )
  }
  
  message("Full DOY total analysis and export complete!")
  return(invisible(NULL))
}

#' Run both standard and total-normalized DOY analyses
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param run_cumulative Whether to run cumulative analysis
#' @param export_data Whether to export data to CSV
#' @return Invisibly returns NULL
run_complete_doy_analysis <- function(years, watersheds, 
                                      run_cumulative = TRUE,
                                      export_data = TRUE) {
  # First run standard DOY analysis
  message("=== Running standard DOY analysis ===")
  run_full_doy_analysis(years, watersheds, run_cumulative, export_data)
  
  # Then run total-normalized DOY analysis
  message("\n=== Running total-normalized DOY analysis ===")
  run_full_doy_total_analysis(years, watersheds, run_cumulative, export_data)
  
  message("\n=== Complete DOY analysis workflow finished ===")
  return(invisible(NULL))
}

# Example usage:

# Run standard DOY timing analysis for Kuskokwim watershed
# run_full_doy_analysis(
#   years = c("2017", "2019", "2020", "2021"),
#   watersheds = c("Kusko"),
#   run_cumulative = TRUE,
#   export_data = TRUE
# )

# Run total-normalized DOY analysis for Kuskokwim watershed
run_full_doy_total_analysis(
  years = c("2017", "2019", "2020", "2021"),
  watersheds = c("Kusko"),
  run_cumulative = TRUE,
  export_data = TRUE
)

# Run both standard and total-normalized analyses together
# run_complete_doy_analysis(
#   years = c("2017", "2019", "2020", "2021"),
#   watersheds = c("Kusko"),
#   run_cumulative = TRUE,
#   export_data = TRUE
# )

# Run DOY timing analysis for both watersheds
# run_complete_doy_analysis(
#   years = c("2015", "2016"),
#   watersheds = c("Yukon", "Kusko"),
#   run_cumulative = TRUE,
#   export_data = TRUE
# )