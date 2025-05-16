# doy_analysis.R
# Functions for DOY (Day of Year) quartile analysis

library(sf)
library(dplyr)
library(here)
library(ggplot2)

# Source the required utility files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

#' Perform DOY Quartile Analysis
#'
#' @param year Character or numeric representing the year
#' @param watershed Character: "Kusko" or "Yukon"
#' @param sensitivity_threshold Numeric threshold for assignment filtering
#' @param min_error Minimum error value to use
#' @param min_stream_order Minimum stream order to include
#' @param HUC HUC level (e.g., 8, 10)
#' @param return_values Whether to return the calculated values
#' @return If return_values is TRUE, a list with results; otherwise NULL
DOY_Quartile_Analysis <- function(year, watershed, sensitivity_threshold, min_error, 
                                  min_stream_order = 3, HUC = 8, 
                                  return_values = FALSE) {
  
  # Generate identifier for output files
  identifier <- paste(year, watershed, sep = "_")
  
  # Load spatial data
  spatial_data <- load_spatial_data(watershed, HUC, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  Huc <- spatial_data$Huc
  
  # Load natal origins data
  natal_data <- load_natal_data(year, watershed)
  
  # Divide data into DOY quartiles
  quartile_data <- divide_doy_quartiles(natal_data)
  quartile_subsets <- quartile_data$subsets
  subset_labels <- quartile_data$labels
  
  # Extract isoscape prediction and error values
  pid_iso <- edges$iso_pred
  pid_isose <- edges$isose_pred
  
  # Calculate error values
  error <- calculate_error(pid_isose, min_error)
  
  # Set up watershed-specific priors
  priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
  
  # Create output directories
  dir.create(here("Basin Maps/DOY_Quartile/HUC"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("Basin Maps/DOY_Quartile/HUC/RawProduction"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("Basin Maps/DOY_Quartile/Tribs"), showWarnings = FALSE, recursive = TRUE)
  
  # Return values storage
  if (return_values) {
    all_results <- list()
  }
  
  # Process each quartile subset
  for (q in 1:length(quartile_subsets)) {
    current_subset <- quartile_subsets[[q]]
    
    # Skip empty subsets
    if (nrow(current_subset) == 0) {
      message(paste("Skipping", subset_labels[q], "because it contains no data"))
      next
    }
    
    # Create unique ID for this subset
    subset_id <- paste0(watershed, "_", year, "_DOY_Q", q)
    message(paste("Processing", subset_labels[q], "with", nrow(current_subset), "data points"))
    
    # Perform assignment
    assignment_matrix <- perform_assignment(
      current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Process assignments to get basin-scale values
    basin_results <- process_assignments(assignment_matrix)
    basin_assign_rescale <- basin_results$rescale
    basin_assign_norm <- basin_results$norm
    
    # Create improved histogram
    gg_hist <- create_doy_histogram(natal_data, current_subset, subset_labels[q])
    
    # Process HUC data
    final_result <- process_huc_data(edges, basin, Huc, basin_assign_rescale, HUC)
    
    # Create HUC map with production per km
    huc_filepath <- file.path(here("Basin Maps/DOY_Quartile/HUC"), 
                              paste0(subset_id, "_HUC", HUC, ".png"))
    
    create_doy_huc_map(
      final_result = final_result,
      basin_assign_norm = basin_assign_norm,
      gg_hist = gg_hist,
      year = year,
      watershed = watershed,
      sensitivity_threshold = sensitivity_threshold,
      min_stream_order = min_stream_order,
      HUC = HUC,
      subset_label = subset_labels[q],
      output_filepath = huc_filepath
    )
    
    # Create HUC map with raw production proportion
    raw_huc_filepath <- file.path(here("Basin Maps/DOY_Quartile/HUC/RawProduction"), 
                                  paste0(subset_id, "_RawProd_HUC", HUC, ".png"))
    
    create_doy_raw_production_map(
      final_result = final_result,
      basin_assign_norm = basin_assign_norm,
      gg_hist = gg_hist,
      year = year,
      watershed = watershed,
      sensitivity_threshold = sensitivity_threshold,
      min_stream_order = min_stream_order,
      HUC = HUC,
      subset_label = subset_labels[q],
      output_filepath = raw_huc_filepath
    )
    
    # Create tributary map
    trib_filepath <- file.path(here("Basin Maps/DOY_Quartile/Tribs"), 
                               paste0(subset_id, ".png"))
    
    create_doy_tributary_map(
      basin = basin,
      edges = edges,
      basin_assign_norm = basin_assign_norm,
      StreamOrderPrior = priors$StreamOrderPrior,
      pid_prior = priors$pid_prior,
      gg_hist = gg_hist,
      year = year,
      watershed = watershed,
      sensitivity_threshold = sensitivity_threshold,
      min_stream_order = min_stream_order,
      min_error = min_error,
      subset_label = subset_labels[q],
      output_filepath = trib_filepath
    )
    
    # Store results if needed
    if (return_values) {
      all_results[[q]] <- list(
        subset = current_subset,
        label = subset_labels[q],
        basin_assign_rescale = basin_assign_rescale,
        basin_assign_norm = basin_assign_norm,
        huc_result = final_result
      )
    }
  }
  
  if (return_values) {
    return(all_results)
  } else {
    return(invisible(NULL))
  }
}