# export_management_analysis.R
# Function to export DOY analysis results with management unit data

library(sf)
library(dplyr)
library(here)
library(readr)
library(tidyr)

# Source the required utility files (make sure these use the modified spatial_utils.R)
source(here("code/utils/spatial_utils.R"))  # Modified version with management units
source(here("code/assignment.R"))
source(here("code/doy_analysis.R"))

#' Export DOY analysis results with management unit breakdown
#'
#' @param years Vector of years to process
#' @param watersheds Vector of watersheds to process (focus on "Kusko" for management units)
#' @param output_dir Directory to save output CSV files
#' @return List with paths to created CSV files
export_management_analysis_results <- function(years, watersheds = "Kusko", 
                                               output_dir = here("Analysis_Results/Management_Units")) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize data frames to store combined results
  all_mgmt_data <- NULL
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
  
  # Helper function to extract tributary data with management info
  extract_trib_mgmt_data <- function(trib_data, basin_assign, metadata) {
    if (is.null(trib_data)) return(NULL)
    
    # Add basin_assign values
    if (!is.null(basin_assign)) {
      trib_data$basin_assign_norm <- basin_assign
    }
    
    # Add metadata
    for (field in names(metadata)) {
      trib_data[[field]] <- metadata[[field]]
    }
    
    # Process as SF object if it is one
    if (inherits(trib_data, "sf")) {
      # Create final data frame with selected columns including management info
      result <- data.frame(
        segment_id = if("reachid" %in% colnames(trib_data)) trib_data$reachid else 1:nrow(trib_data),
        stream_order = if("Str_Order" %in% colnames(trib_data)) trib_data$Str_Order else NA,
        mgmt_river = if("mgmt_river" %in% colnames(trib_data)) trib_data$mgmt_river else NA,
        stream_length_m = as.numeric(st_length(trib_data)),
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
      message(paste("Processing management unit analysis for", watershed, "watershed, year", year))
      
      # Set parameters based on watershed
      if (watershed == "Yukon") {
        sensitivity_threshold <- 0.7
        min_error <- 0.003
        min_stream_order <- 5
        message("Note: Management unit analysis is primarily designed for Kusko watershed")
      } else if (watershed == "Kusko") {
        sensitivity_threshold <- 0.7
        min_error <- 0.0006
        min_stream_order <- 3
      } else {
        warning(paste("Unknown watershed:", watershed, "- skipping"))
        next
      }
      
      # Load spatial data (using modified function that loads management shapefile)
      spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
      edges <- spatial_data$edges
      basin <- spatial_data$basin
      Huc <- spatial_data$Huc
      
      # Check if we have management unit data
      if (!"mgmt_river" %in% colnames(edges)) {
        warning(paste("mgmt_river column not found for", watershed, "- skipping management analysis"))
        next
      }
      
      # Load natal origins data
      natal_data <- load_natal_data(year, watershed)
      
      # Divide data into DOY quartiles
      quartile_data <- divide_doy_quartiles(natal_data)
      quartile_subsets <- quartile_data$subsets
      
      # Extract isoscape prediction and error values
      pid_iso <- edges$iso_pred
      pid_isose <- edges$isose_pred
      
      # Calculate error values
      error <- calculate_error(pid_isose, min_error)
      
      # Set up watershed-specific priors
      priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
      
      # First, calculate TOTAL production across all quartiles for normalization
      message("  Calculating total production across all quartiles...")
      all_data <- do.call(rbind, quartile_subsets)
      all_assignment_matrix <- perform_assignment(
        all_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
      )
      
      # Calculate total basin production
      total_basin_assign_sum <- apply(all_assignment_matrix, 1, sum, na.rm = TRUE)
      grand_total_production <- sum(total_basin_assign_sum, na.rm = TRUE)
      
      # Calculate total management unit proportions for the entire run
      total_mgmt_result <- process_management_data(edges, total_basin_assign_sum, grand_total_production)
      if (is.null(total_mgmt_result)) {
        warning(paste("No management unit data available for", watershed, "in", year))
        next
      }
      
      message(paste("  Total production:", round(grand_total_production, 6)))
      managed_total <- sum(total_mgmt_result$total_production, na.rm = TRUE)
      message(paste("  Managed production:", round(managed_total, 6), 
                    "| Managed % of total:", round(100 * managed_total / grand_total_production, 1), "%"))
      
      # Process each DOY quartile
      for (q in 1:length(quartile_subsets)) {
        current_subset <- quartile_subsets[[q]]
        
        if (nrow(current_subset) == 0) {
          message(paste("    Skipping Q", q, "because it contains no data"))
          next
        }
        
        message(paste("    Processing DOY Q", q, "with", nrow(current_subset), "data points"))
        
        # Perform assignment for this quartile
        assignment_matrix <- perform_assignment(
          current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
        )
        
        # Calculate basin values for this quartile
        quartile_basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
        
        # Process management unit data for this quartile
        mgmt_result <- process_management_data(edges, quartile_basin_assign_sum, grand_total_production)
        
        # Check if we got any results
        if (is.null(mgmt_result) || nrow(mgmt_result) == 0) {
          message(paste("    Warning: No management unit results for Q", q, "- skipping"))
          next
        }
        
        # mgmt_result already has percent_of_total_run calculated correctly
        # Add some additional useful metrics
        quartile_managed_total <- sum(mgmt_result$total_production, na.rm = TRUE)
        mgmt_result$quartile_managed_total <- quartile_managed_total
        mgmt_result$quartile_managed_percent_of_total <- (quartile_managed_total / grand_total_production) * 100
        
        # Add metadata
        mgmt_metadata <- list(
          year = year,
          watershed = watershed,
          analysis_type = "DOY_Management",
          quartile = paste0("Q", q)
        )
        
        # Add management unit data
        mgmt_df <- mgmt_result
        for (field in names(mgmt_metadata)) {
          mgmt_df[[field]] <- mgmt_metadata[[field]]
        }
        
        # Debug: Print info about mgmt_df
        message(paste("    Management units found:", nrow(mgmt_df)))
        if (nrow(mgmt_df) > 0) {
          message(paste("    Management rivers:", paste(unique(mgmt_df$mgmt_river), collapse = ", ")))
          message(paste("    Quartile managed % of total run:", 
                        round(quartile_managed_total / grand_total_production * 100, 2), "%"))
        }
        
        all_mgmt_data <- rbind(all_mgmt_data, mgmt_df)
        
        # Also process HUC data for comparison
        huc_result <- process_huc_data(edges, basin, Huc, quartile_basin_assign_sum, 8)
        huc_result$percent_of_total_run <- (huc_result$total_production / grand_total_production) * 100
        
        huc_metadata <- list(
          year = year,
          watershed = watershed,
          analysis_type = "DOY_HUC",
          quartile = paste0("Q", q)
        )
        
        huc_df <- extract_data(huc_result)
        for (field in names(huc_metadata)) {
          huc_df[[field]] <- huc_metadata[[field]]
        }
        all_huc_data <- rbind(all_huc_data, huc_df)
        
        # Add tributary data with management info
        trib_metadata <- list(
          year = year,
          watershed = watershed,
          analysis_type = "DOY_Tributary",
          quartile = paste0("Q", q)
        )
        
        # Calculate normalized assignments for tributaries
        quartile_basin_assign_norm <- quartile_basin_assign_sum / max(quartile_basin_assign_sum, na.rm = TRUE)
        
        trib_df <- extract_trib_mgmt_data(
          edges, 
          quartile_basin_assign_norm, 
          trib_metadata
        )
        all_trib_data <- rbind(all_trib_data, trib_df)
      }
    }
  }
  
  # Save management unit data
  if (!is.null(all_mgmt_data) && nrow(all_mgmt_data) > 0) {
    mgmt_filepath <- file.path(output_dir, "management_unit_production_by_quartile.csv")
    write.csv(all_mgmt_data, mgmt_filepath, row.names = FALSE)
    message(paste("Saved management unit data to:", mgmt_filepath))
  } else {
    warning("No management unit data collected")
    mgmt_filepath <- NULL
  }
  
  # Save HUC data
  if (!is.null(all_huc_data) && nrow(all_huc_data) > 0) {
    huc_filepath <- file.path(output_dir, "huc_production_by_quartile_with_mgmt.csv")
    write.csv(all_huc_data, huc_filepath, row.names = FALSE)
    message(paste("Saved HUC data to:", huc_filepath))
  } else {
    warning("No HUC data collected")
    huc_filepath <- NULL
  }
  
  # Save tributary data with management info
  if (!is.null(all_trib_data) && nrow(all_trib_data) > 0) {
    trib_filepath <- file.path(output_dir, "tributary_production_by_quartile_with_mgmt.csv")
    write.csv(all_trib_data, trib_filepath, row.names = FALSE)
    message(paste("Saved tributary data to:", trib_filepath))
  } else {
    warning("No tributary data collected")
    trib_filepath <- NULL
  }
  
  # Create a summary table showing management unit contributions across quartiles
  if (!is.null(all_mgmt_data) && nrow(all_mgmt_data) > 0) {
    summary_table <- all_mgmt_data %>%
      select(year, quartile, mgmt_river, 
             percent_of_total_run, total_production, edge_count) %>%
      arrange(year, quartile, desc(percent_of_total_run))
    
    summary_filepath <- file.path(output_dir, "management_unit_summary.csv")
    write.csv(summary_table, summary_filepath, row.names = FALSE)
    message(paste("Saved management unit summary to:", summary_filepath))
  } else {
    summary_filepath <- NULL
  }
  
  return(list(
    mgmt_file = mgmt_filepath,
    huc_file = huc_filepath,
    trib_file = trib_filepath,
    summary_file = summary_filepath
  ))
}

#' Create a visualization of management unit contributions by quartile
#'
#' @param mgmt_data_file Path to the management unit CSV file
#' @param output_dir Directory to save the plot
#' @return Path to the created plot
create_management_unit_plot <- function(mgmt_data_file, output_dir = here("Analysis_Results/Management_Units")) {
  
  if (!file.exists(mgmt_data_file)) {
    warning("Management unit data file not found")
    return(NULL)
  }
  
  # Load the data
  mgmt_data <- read.csv(mgmt_data_file)
  
  # Create a plot showing management unit contributions by quartile and year
  library(ggplot2)
  library(RColorBrewer)
  
  # Filter out any management units with very small contributions for clarity
  mgmt_data_filtered <- mgmt_data %>%
    group_by(mgmt_river) %>%
    filter(max(percent_of_total_run) >= 0.1) %>%  # Only include if max contribution >= 0.1%
    ungroup()
  
  # Create the plot
  p <- ggplot(mgmt_data_filtered, aes(x = quartile, y = percent_of_total_run, fill = mgmt_river)) +
    geom_col(position = "dodge", alpha = 0.8) +
    facet_wrap(~year, scales = "free_x") +
    scale_fill_brewer(palette = "Set3", name = "Management\nRiver") +
    scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
    labs(
      title = "Management Unit Contributions to Total Salmon Run by DOY Quartile",
      subtitle = "Percentage of ENTIRE watershed run by management river and timing quartile",
      x = "DOY Quartile",
      y = "Percent of Total Watershed Run",
      caption = "Note: Many tributaries are unmanaged and not shown. Only management units contributing ≥0.1% shown."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  # Save the plot
  plot_path <- file.path(output_dir, "management_unit_contributions_by_quartile.png")
  ggsave(plot_path, p, width = 12, height = 8, dpi = 300)
  
  message(paste("Created management unit plot:", plot_path))
  return(plot_path)
}

# Example usage function
#' Run the complete management unit analysis
#'
#' @param years Vector of years to analyze
#' @param create_plot Whether to create visualization
#' @return List of output file paths
run_management_analysis <- function(years = c("2017", "2019", "2020", "2021"), 
                                    create_plot = TRUE) {
  
  message("Starting management unit analysis...")
  
  # Export the analysis results
  results <- export_management_analysis_results(
    years = years,
    watersheds = "Kusko"  # Focus on Kusko since that's where we have management units
  )
  
  # Create visualization if requested and we have data
  if (create_plot && !is.null(results$mgmt_file)) {
    plot_path <- create_management_unit_plot(results$mgmt_file)
    results$plot_file <- plot_path
  }
  
  message("Management unit analysis complete!")
  return(results)
}

# Uncomment to run the analysis:
 run_management_analysis(years = c("2017", "2019", "2020", "2021"))
 