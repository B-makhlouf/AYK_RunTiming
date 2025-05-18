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

#' Export DOY coefficient of variation analysis results
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process
#' @param quartiles Vector of quartiles to analyze (e.g., c("Q1", "Q2", "Q3", "Q4"))
#' @param top_percent Percent of top producers to include in CV analysis
#' @param output_dir Directory to save output CSV files
#' @return List with paths to created CSV files
export_doy_cv_analysis <- function(years, watersheds, 
                                   quartiles = c("Q1", "Q2", "Q3", "Q4"),
                                   top_percent = 50,
                                   output_dir = here("Analysis_Results")) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize data frames to store combined results
  all_cv_data <- NULL
  
  # Process each watershed
  for (watershed in watersheds) {
    message(paste("Processing CV for", watershed, "watershed"))
    
    # Set watershed-specific parameters
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
    
    # Process each quartile
    for (quartile in quartiles) {
      message(paste("  Processing", quartile))
      
      tryCatch({
        # Collect data across years for this quartile
        quartile_data_list <- list()
        
        for (year in years) {
          # Run DOY analysis for this year
          doy_results <- DOY_Quartile_Analysis(
            year = year,
            watershed = watershed,
            sensitivity_threshold = sensitivity_threshold,
            min_error = min_error,
            min_stream_order = min_stream_order,
            HUC = 8,
            return_values = TRUE
          )
          
          # Extract the specific quartile data
          q_index <- as.numeric(gsub("Q", "", quartile))
          if (q_index <= length(doy_results) && !is.null(doy_results[[q_index]])) {
            quartile_data_list[[year]] <- list(
              huc_result = doy_results[[q_index]]$huc_result,
              basin_assign_norm = doy_results[[q_index]]$basin_assign_norm
            )
          }
        }
        
        # Calculate CV across years for this quartile
        if (length(quartile_data_list) > 1) {
          cv_result <- calculate_doy_cv_across_years(
            quartile_data_list, 
            watershed,
            quartile,
            top_percent
          )
          
          # Add metadata
          cv_result$cv_data$watershed <- watershed
          cv_result$cv_data$analysis_type <- paste0("DOY_", quartile)
          cv_result$cv_data$top_percent <- top_percent
          
          # Combine with all CV data
          all_cv_data <- rbind(all_cv_data, cv_result$cv_data)
        }
        
      }, error = function(e) {
        message(paste("    Error processing", quartile, "CV:", e$message))
      })
    }
  }
  
  # Save CV data to CSV
  if (!is.null(all_cv_data) && nrow(all_cv_data) > 0) {
    cv_filepath <- file.path(output_dir, paste0("doy_cv_analysis_top", top_percent, ".csv"))
    write.csv(all_cv_data, cv_filepath, row.names = FALSE)
    message(paste("Saved DOY CV analysis to:", cv_filepath))
    
    return(list(cv_file = cv_filepath))
  } else {
    warning("No CV data collected - unable to save CV CSV")
    return(list(cv_file = NULL))
  }
}

#' Calculate coefficient of variation across years for DOY analysis
#'
#' @param quartile_data_list List of yearly data for a specific quartile
#' @param watershed Character: "Kusko" or "Yukon"
#' @param quartile Character: e.g., "Q1", "Q2", "Q3", "Q4"
#' @param top_percent Percent of top producers to include
#' @return List with CV analysis results
calculate_doy_cv_across_years <- function(quartile_data_list, watershed, quartile, top_percent) {
  # Extract HUC data from each year
  years <- names(quartile_data_list)
  n_years <- length(years)
  
  # Get reference HUC structure from first year
  ref_huc <- quartile_data_list[[1]]$huc_result
  ref_huc <- st_drop_geometry(ref_huc)
  
  n_hucs <- nrow(ref_huc)
  
  # Initialize matrices to store values across years
  production_matrix <- matrix(NA, nrow = n_hucs, ncol = n_years)
  
  # Fill matrix with values from each year
  for (i in 1:n_years) {
    year <- years[i]
    current_huc <- quartile_data_list[[year]]$huc_result
    current_huc <- st_drop_geometry(current_huc)
    
    # Match HUCs between reference and current year
    match_idx <- match(ref_huc$HUC8, current_huc$HUC8)
    
    # Extract production values
    production_matrix[, i] <- current_huc$production_proportion[match_idx]
  }
  
  # Calculate mean production across years
  mean_production <- rowMeans(production_matrix, na.rm = TRUE)
  
  # Identify top producers based on mean production
  production_threshold <- quantile(mean_production, prob = 1 - (top_percent/100), na.rm = TRUE)
  top_producers <- mean_production >= production_threshold
  
  # Calculate CV only for top producers
  cv_values <- rep(NA, n_hucs)
  
  for (i in 1:n_hucs) {
    if (top_producers[i] && sum(!is.na(production_matrix[i, ])) >= 2) {
      mean_val <- mean_production[i]
      if (mean_val > 0) {
        cv_values[i] <- sd(production_matrix[i, ], na.rm = TRUE) / mean_val
      }
    }
  }
  
  # Create result data frame
  cv_data <- data.frame(
    HUC8 = ref_huc$HUC8,
    Name = ref_huc$Name,
    mean_production = mean_production,
    cv_production = cv_values,
    is_top_producer = top_producers,
    n_years = rowSums(!is.na(production_matrix)),
    quartile = quartile
  )
  
  return(list(
    cv_data = cv_data,
    production_matrix = production_matrix,
    years = years
  ))
}