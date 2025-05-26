# doy_analysis.R - FIXED VERSION
# Functions for DOY (Day of Year) quartile analysis

library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(grid)
library(gridExtra)

# Source the required utility files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

#' Perform DOY Quartile Analysis - FIXED VERSION
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
  
  message(paste("=== Starting DOY Quartile Analysis for", year, watershed, "==="))
  
  # Generate identifier for output files
  identifier <- paste(year, watershed, sep = "_")
  
  # Load spatial data
  message("Loading spatial data...")
  tryCatch({
    spatial_data <- load_spatial_data(watershed, HUC, min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    Huc <- spatial_data$Huc
    message(paste("  - Loaded", nrow(edges), "edges,", nrow(basin), "basin features,", nrow(Huc), "HUC features"))
  }, error = function(e) {
    stop("Error loading spatial data: ", e$message)
  })
  
  # Load natal origins data
  message("Loading natal origins data...")
  tryCatch({
    natal_data <- load_natal_data(year, watershed)
    message(paste("  - Loaded", nrow(natal_data), "natal data records"))
  }, error = function(e) {
    stop("Error loading natal data: ", e$message)
  })
  
  # Divide data into DOY quartiles
  message("Dividing data into DOY quartiles...")
  tryCatch({
    quartile_data <- divide_doy_quartiles(natal_data)
    quartile_subsets <- quartile_data$subsets
    subset_labels <- quartile_data$labels
    
    for (q in 1:length(quartile_subsets)) {
      message(paste("  - Quartile", q, ":", nrow(quartile_subsets[[q]]), "records"))
    }
  }, error = function(e) {
    stop("Error dividing into quartiles: ", e$message)
  })
  
  # Extract isoscape prediction and error values
  message("Setting up predictions and error calculations...")
  tryCatch({
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    
    # Calculate error values
    error <- calculate_error(pid_isose, min_error)
    
    # Set up watershed-specific priors
    priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
    
    message(paste("  - Setup complete with", length(pid_iso), "prediction points"))
  }, error = function(e) {
    stop("Error setting up predictions: ", e$message)
  })
  
  # Create output directories
  message("Creating output directories...")
  dir.create(here("Basin Maps/DOY_Quartile/HUC"), showWarnings = FALSE, recursive = TRUE)
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
      message(paste("  Skipping", subset_labels[q], "- no data"))
      next
    }
    
    message(paste("Processing", subset_labels[q], "with", nrow(current_subset), "data points..."))
    
    # Create unique ID for this subset
    subset_id <- paste0(watershed, "_", year, "_DOY_Q", q)
    
    # Perform assignment
    tryCatch({
      assignment_matrix <- perform_assignment(
        current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
      )
      
      # Process assignments to get basin-scale values
      basin_results <- process_assignments(assignment_matrix)
      basin_assign_rescale <- basin_results$rescale
      basin_assign_norm <- basin_results$norm
      
      message(paste("  - Assignment completed, max value:", round(max(basin_assign_norm, na.rm = TRUE), 4)))
      
    }, error = function(e) {
      message(paste("  Error in assignment for quartile", q, ":", e$message))
      next
    })
    
    # Create improved histogram
    tryCatch({
      gg_hist <- create_doy_histogram(natal_data, current_subset, subset_labels[q])
    }, error = function(e) {
      message(paste("  Warning: Could not create histogram:", e$message))
      gg_hist <- NULL
    })
    
    # Process HUC data
    tryCatch({
      final_result <- process_huc_data(edges, basin, Huc, basin_assign_rescale, HUC)
      message(paste("  - HUC processing completed for", nrow(final_result), "HUCs"))
    }, error = function(e) {
      message(paste("  Error processing HUC data for quartile", q, ":", e$message))
      next
    })
    
    # Create HUC map
    huc_filepath <- file.path(here("Basin Maps/DOY_Quartile/HUC"), 
                              paste0(subset_id, "_HUC", HUC, ".png"))
    
    message(paste("  - Creating HUC map:", basename(huc_filepath)))
    tryCatch({
      # Create HUC map using the visualization function
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
      
      if (file.exists(huc_filepath)) {
        message(paste("    ✓ HUC map created successfully"))
      } else {
        message(paste("    ✗ HUC map file not found after creation"))
      }
      
    }, error = function(e) {
      message(paste("  Error creating HUC map for quartile", q, ":", e$message))
    })
    
    # Create tributary map
    trib_filepath <- file.path(here("Basin Maps/DOY_Quartile/Tribs"), 
                               paste0(subset_id, ".png"))
    
    message(paste("  - Creating tributary map:", basename(trib_filepath)))
    tryCatch({
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
      
      if (file.exists(trib_filepath)) {
        message(paste("    ✓ Tributary map created successfully"))
      } else {
        message(paste("    ✗ Tributary map file not found after creation"))
      }
      
    }, error = function(e) {
      message(paste("  Error creating tributary map for quartile", q, ":", e$message))
    })
    
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
    
    message(paste("  ✓ Completed processing quartile", q))
  }
  
  message(paste("=== DOY Quartile Analysis completed for", year, watershed, "==="))
  
  if (return_values) {
    return(all_results)
  } else {
    return(invisible(NULL))
  }
}

#' Test function to debug DOY analysis issues
#'
#' @param year Character or numeric representing the year
#' @param watershed Character: "Kusko" or "Yukon"
test_doy_quartile_analysis <- function(year = "2019", watershed = "Kusko") {
  message("=== Testing DOY Quartile Analysis ===")
  
  # Set parameters
  if (watershed == "Kusko") {
    sensitivity_threshold <- 0.7
    min_error <- 0.0006
    min_stream_order <- 3
  } else {
    sensitivity_threshold <- 0.7
    min_error <- 0.003
    min_stream_order <- 5
  }
  
  # Test the analysis
  result <- DOY_Quartile_Analysis(
    year = year,
    watershed = watershed,
    sensitivity_threshold = sensitivity_threshold,
    min_error = min_error,
    min_stream_order = min_stream_order,
    HUC = 8,
    return_values = FALSE
  )
  
  # Check if files were created
  output_dir_huc <- here("Basin Maps/DOY_Quartile/HUC")
  output_dir_tribs <- here("Basin Maps/DOY_Quartile/Tribs")
  
  huc_files <- list.files(output_dir_huc, pattern = paste0(watershed, "_", year, "_DOY"), full.names = TRUE)
  trib_files <- list.files(output_dir_tribs, pattern = paste0(watershed, "_", year, "_DOY"), full.names = TRUE)
  
  message(paste("Files created:"))
  message(paste("  HUC maps:", length(huc_files)))
  message(paste("  Tributary maps:", length(trib_files)))
  
  if (length(huc_files) > 0) {
    message("  HUC files:")
    for (f in huc_files) message(paste("   ", basename(f)))
  }
  
  if (length(trib_files) > 0) {
    message("  Tributary files:")
    for (f in trib_files) message(paste("   ", basename(f)))
  }
  
  return(invisible(NULL))
}

# Example usage to test:
# test_doy_quartile_analysis("2019", "Kusko")