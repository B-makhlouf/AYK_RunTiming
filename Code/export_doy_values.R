# export_doy_values.R
# Function to export DOY analysis results to CSV files

library(sf)
library(dplyr)
library(here)
library(readr)

#' Export DOY analysis results to CSV files
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param include_cumulative Whether to include cumulative DOY analysis
#' @param output_dir Directory to save output CSV files
#' @param cache_dir Directory to look for cached analysis results
#' @return List with paths to created CSV files
export_doy_analysis_results <- function(years, watersheds, 
                                        include_cumulative = TRUE,
                                        output_dir = here("Analysis_Results"),
                                        cache_dir = here("Data/Processed")) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize data frames to store combined results
  all_huc_data <- NULL
  all_trib_data <- NULL
  
  # Helper function to extract non-geometry data from SF objects
  extract_data <- function(sf_obj) {
    if (inherits(sf_obj, "sf")) {
      return(st_drop_geometry(sf_obj))
    } else {
      return(sf_obj)
    }
  }
  
  # Helper function to extract tributary data with consistent columns
  extract_trib_data <- function(trib_data, basin_assign, metadata) {
    # Early return if input is NULL
    if (is.null(trib_data)) return(NULL)
    
    # Add basin_assign_norm if provided
    if (!is.null(basin_assign)) {
      trib_data$basin_assign_norm <- basin_assign
    }
    
    # Add metadata
    for (field in names(metadata)) {
      trib_data[[field]] <- metadata[[field]]
    }
    
    # Process as SF object if it is one
    if (inherits(trib_data, "sf")) {
      # Determine ID column
      id_col <- NULL
      for (possible_id in c("reachid", "SEGMENT_ID", "REACH_ID", "id")) {
        if (possible_id %in% colnames(trib_data)) {
          id_col <- possible_id
          break
        }
      }
      
      # Create ID if none exists
      if (is.null(id_col)) {
        trib_data$segment_id <- 1:nrow(trib_data)
        id_col <- "segment_id"
      } else {
        # Rename to standard column name
        trib_data$segment_id <- trib_data[[id_col]]
      }
      
      # Get stream order
      order_col <- NULL
      for (possible_order in c("Str_Order", "STR_ORDER", "stream_order")) {
        if (possible_order %in% colnames(trib_data)) {
          order_col <- possible_order
          break
        }
      }
      
      # Default stream order if none exists
      if (is.null(order_col)) {
        trib_data$stream_order <- NA
      } else {
        # Rename to standard column name
        trib_data$stream_order <- trib_data[[order_col]]
      }
      
      # Add stream length if not present
      if (!"stream_length_m" %in% colnames(trib_data)) {
        trib_data$stream_length_m <- as.numeric(st_length(trib_data))
      }
      
      # Create final data frame with selected columns
      result <- data.frame(
        segment_id = trib_data$segment_id,
        stream_order = trib_data$stream_order,
        stream_length_m = trib_data$stream_length_m,
        basin_assign_norm = trib_data$basin_assign_norm
      )
      
      # Add metadata columns
      for (field in names(metadata)) {
        result[[field]] <- metadata[[field]]
      }
      
      return(result)
    } else {
      # Already a data frame, return as is
      return(trib_data)
    }
  }
  
  # Process each watershed and year
  for (watershed in watersheds) {
    for (year in years) {
      message(paste("Processing DOY analysis results for", watershed, "watershed, year", year))
      
      # Set parameters based on watershed
      if (watershed == "Yukon") {
        sensitivity_threshold <- 0.7
        min_error <- 0.003
        min_stream_order <- 5
      } else if (watershed == "Kusko") {
        sensitivity_threshold <- 0.7
        min_error <- 0.0006
        min_stream_order <- 3
      } else {
        warning(paste("Unknown watershed:", watershed, "- skipping"))
        next
      }
      
      # Load spatial data (needed for all analyses)
      spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
      edges <- spatial_data$edges
      basin <- spatial_data$basin
      Huc <- spatial_data$Huc
      
      # Process DOY quartile data
      message("  Processing DOY quartile data")
      
      tryCatch({
        doy_results <- DOY_Quartile_Analysis(
          year = year,
          watershed = watershed,
          sensitivity_threshold = sensitivity_threshold,
          min_error = min_error,
          min_stream_order = min_stream_order,
          HUC = 8,
          return_values = TRUE
        )
        
        # Process each DOY quartile
        for (q in 1:length(doy_results)) {
          if (is.null(doy_results[[q]])) next
          
          message(paste("    Processing DOY Q", q))
          
          doy_metadata <- list(
            year = year,
            watershed = watershed,
            analysis_type = "DOY",
            quartile = paste0("Q", q)
          )
          
          # Add HUC data
          doy_huc_df <- extract_data(doy_results[[q]]$huc_result)
          for (field in names(doy_metadata)) {
            doy_huc_df[[field]] <- doy_metadata[[field]]
          }
          all_huc_data <- rbind(all_huc_data, doy_huc_df)
          
          # Add tributary data
          doy_trib_df <- extract_trib_data(
            edges, 
            doy_results[[q]]$basin_assign_norm, 
            doy_metadata
          )
          all_trib_data <- rbind(all_trib_data, doy_trib_df)
        }
      }, error = function(e) {
        message("    Error processing DOY quartiles: ", e$message)
      })
      
      # Process cumulative DOY data if requested
      if (include_cumulative) {
        message("  Processing cumulative DOY data")
        
        tryCatch({
          cum_doy_results <- Cumulative_DOY_Analysis(
            year = year,
            watershed = watershed,
            sensitivity_threshold = sensitivity_threshold,
            min_error = min_error,
            min_stream_order = min_stream_order,
            HUC = 8,
            return_values = TRUE
          )
          
          # Process each cumulative DOY quartile
          for (q in 1:length(cum_doy_results)) {
            if (is.null(cum_doy_results[[q]])) next
            
            message(paste("    Processing Cumulative DOY Q", q))
            
            cum_doy_metadata <- list(
              year = year,
              watershed = watershed,
              analysis_type = "Cumulative_DOY",
              quartile = paste0("Q", q)
            )
            
            # Add HUC data
            cum_doy_huc_df <- extract_data(cum_doy_results[[q]]$huc_result)
            for (field in names(cum_doy_metadata)) {
              cum_doy_huc_df[[field]] <- cum_doy_metadata[[field]]
            }
            all_huc_data <- rbind(all_huc_data, cum_doy_huc_df)
            
            # Add tributary data
            cum_doy_trib_df <- extract_trib_data(
              edges, 
              cum_doy_results[[q]]$basin_assign_norm, 
              cum_doy_metadata
            )
            all_trib_data <- rbind(all_trib_data, cum_doy_trib_df)
          }
        }, error = function(e) {
          message("    Error processing cumulative DOY data: ", e$message)
        })
      }
    }
  }
  
  # Check if we have data to save
  if (is.null(all_huc_data) || nrow(all_huc_data) == 0) {
    warning("No HUC data collected - unable to save HUC CSV")
    huc_filepath <- NULL
  } else {
    # Save HUC data to CSV
    huc_filepath <- file.path(output_dir, "all_doy_huc_values.csv")
    write.csv(all_huc_data, huc_filepath, row.names = FALSE)
    message(paste("Saved HUC data to:", huc_filepath))
  }
  
  if (is.null(all_trib_data) || nrow(all_trib_data) == 0) {
    warning("No tributary data collected - unable to save tributary CSV")
    trib_filepath <- NULL
  } else {
    # Save tributary data to CSV
    trib_filepath <- file.path(output_dir, "all_doy_tributary_values.csv")
    write.csv(all_trib_data, trib_filepath, row.names = FALSE)
    message(paste("Saved tributary data to:", trib_filepath))
  }
  
  return(list(
    huc_file = huc_filepath,
    trib_file = trib_filepath
  ))
}


export_doy_analysis_results(
  watersheds = "Kusko",
  years = 2018,
  include_cumulative = FALSE,
  output_dir = "output"
)
