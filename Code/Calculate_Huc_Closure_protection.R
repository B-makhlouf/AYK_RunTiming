# huc_closure_proportions.R
# Script to calculate and export HUC proportions during the closure window
# Exports the green vs. gray proportions as seen in the visualization

library(dplyr)
library(sf)
library(here)
library(readr)
library(tidyr) # For data manipulation functions

# Source the necessary utility files
source(here("code/utils/spatial_utils.R"))
source(here("code/assignment.R"))

#' Calculate and export HUC proportions during closure window
#'
#' @param years Vector of years to analyze
#' @param watershed Character: "Kusko" or "Yukon" 
#' @param closure_start_date Character: Closure window start date in "MM-DD" format (default "06-01")
#' @param closure_end_date Character: Closure window end date in "MM-DD" format (default "06-11")
#' @param output_dir Directory to save output CSV files
#' @param HUC HUC level (e.g., 8, 10)
#' @return Path to created CSV file
export_huc_closure_proportions <- function(years, 
                                           watershed, 
                                           closure_start_date = "06-01", 
                                           closure_end_date = "06-11",
                                           output_dir = here("Analysis_Results"),
                                           HUC = 8) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize data frame to store results
  all_results <- data.frame()
  
  # Convert closure dates to DOY (Day of Year)
  # Note: Use a leap year to ensure proper DOY calculation
  closure_start_doy <- as.numeric(format(as.Date(paste0("2024-", closure_start_date)), "%j"))
  closure_end_doy <- as.numeric(format(as.Date(paste0("2024-", closure_end_date)), "%j"))
  
  message(paste("Analyzing closure window DOY", closure_start_doy, "to", closure_end_doy, 
                "(", closure_start_date, "to", closure_end_date, ")"))
  
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
    stop(paste("Unknown watershed:", watershed))
  }
  
  # Process each year
  for (year in years) {
    tryCatch({
      message(paste("Processing", watershed, "watershed for year", year))
      
      # Load spatial data
      spatial_data <- load_spatial_data(watershed, HUC, min_stream_order)
      edges <- spatial_data$edges
      basin <- spatial_data$basin
      Huc <- spatial_data$Huc
      
      # Load natal origins data
      natal_data <- load_natal_data(year, watershed)
      
      # Filter the natal data to only include fish within the closure window
      window_data <- natal_data %>% 
        filter(DOY >= closure_start_doy, DOY <= closure_end_doy)
      
      if (nrow(window_data) == 0) {
        warning(paste("No data found in closure window for", watershed, "in", year))
        next
      }
      
      message(paste("  Found", nrow(window_data), "of", nrow(natal_data), 
                    "total observations within closure window"))
      
      # Extract isoscape prediction and error values
      pid_iso <- edges$iso_pred
      pid_isose <- edges$isose_pred
      
      # Calculate error values
      error <- calculate_error(pid_isose, min_error)
      
      # Set up watershed-specific priors
      priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
      
      # --- Calculate total production across all DOY values ---
      message("  Calculating total production across all DOY values...")
      
      total_assignment_matrix <- perform_assignment(
        natal_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
      )
      
      # Process total assignments to get basin-scale values
      total_basin_results <- process_assignments(total_assignment_matrix)
      total_basin_assign_sum <- total_basin_results$sum
      
      # --- Create a custom version of process_huc_data that doesn't use replace_na ---
      custom_process_huc_data <- function(edges, basin, Huc, basin_assign_rescale, HUC = 8) {
        huc_col <- paste0("HUC", HUC)
        name_col <- "Name"
        
        # remove the column HUC8 from edges 
        if (huc_col %in% colnames(edges)) {
          edges <- edges %>% select(-all_of(huc_col))
        }
        
        # Transform to consistent CRS
        edges <- st_transform(edges, st_crs(Huc))
        edges$basin_assign_rescale <- basin_assign_rescale
        basin <- st_transform(basin, st_crs(Huc))
        
        # Identify which HUCs intersect with the basin
        basin_buffer <- st_buffer(basin, dist = 0)
        hucs_in_basin <- Huc[st_intersects(Huc, basin_buffer, sparse = FALSE)[,1], ]
        
        # Calculate intersection areas
        intersection_areas <- st_intersection(hucs_in_basin, basin_buffer) %>%
          mutate(area = st_area(.)) %>%
          st_drop_geometry() %>%
          group_by(!!sym(huc_col)) %>%
          summarize(int_area = sum(area))
        
        # Get original HUC areas
        hucs_areas <- hucs_in_basin %>%
          mutate(total_area = st_area(.)) %>%
          st_drop_geometry() %>%
          select(!!sym(huc_col), total_area)
        
        # Calculate overlap percentage
        overlap_percentage <- intersection_areas %>%
          left_join(hucs_areas, by = huc_col) %>%
          mutate(pct_overlap = as.numeric(int_area / total_area))
        
        # Filter to include only HUCs with significant overlap
        significant_hucs <- overlap_percentage %>%
          filter(pct_overlap > 0.1) %>%
          pull(!!sym(huc_col))
        
        # Spatial join between edges and HUCs
        Combined_edges_HUC <- st_join(edges, Huc, join = st_intersects)
        edges$stream_length_m <- as.numeric(st_length(edges))
        
        # Summarize production by HUC polygon
        summary_huc <- Combined_edges_HUC %>%
          group_by(!!sym(huc_col)) %>%
          summarise(
            total_production = sum(basin_assign_rescale, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Replace NAs with 0 in production values
        summary_huc$total_production[is.na(summary_huc$total_production)] <- 0
        
        # Calculate production proportion (handle case where all are 0)
        total_prod <- sum(summary_huc$total_production, na.rm = TRUE)
        if (total_prod > 0) {
          summary_huc$production_proportion <- summary_huc$total_production / total_prod
        } else {
          summary_huc$production_proportion <- 0
        }
        
        # Summarize stream length by HUC
        stream_length_by_huc <- edges %>%
          st_join(Huc, join = st_intersects) %>%
          st_drop_geometry() %>%
          group_by(!!sym(huc_col)) %>%
          summarise(total_stream_length = sum(stream_length_m, na.rm = TRUE))
        
        # Replace NAs with 0 in stream length
        stream_length_by_huc$total_stream_length[is.na(stream_length_by_huc$total_stream_length)] <- 0
        
        # Merge production and stream length data with HUC polygons
        final_result <- Huc %>%
          filter(!!sym(huc_col) %in% significant_hucs) %>%
          left_join(st_drop_geometry(summary_huc), by = huc_col) %>%
          left_join(stream_length_by_huc, by = huc_col)
        
        # Replace NAs with 0
        final_result$total_production[is.na(final_result$total_production)] <- 0
        final_result$production_proportion[is.na(final_result$production_proportion)] <- 0
        final_result$total_stream_length[is.na(final_result$total_stream_length)] <- 0
        
        # Calculate additional metrics
        final_result <- final_result %>% mutate(
          production_per_meter = ifelse(total_stream_length > 0,
                                        total_production / total_stream_length,
                                        0)
        )
        
        # Calculate normalized production per meter
        max_prod_per_meter <- max(final_result$production_per_meter, na.rm = TRUE)
        if (max_prod_per_meter > 0) {
          final_result$production_per_meter_norm <- final_result$production_per_meter / max_prod_per_meter
        } else {
          final_result$production_per_meter_norm <- 0
        }
        
        return(final_result)
      }
      
      # Process HUC data for total production using custom function
      total_huc_result <- custom_process_huc_data(edges, basin, Huc, total_basin_assign_sum, HUC)
      
      # --- Calculate production within closure window ---
      message("  Calculating production within closure window...")
      
      window_assignment_matrix <- perform_assignment(
        window_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
      )
      
      # Process window assignments to get basin-scale values
      window_basin_results <- process_assignments(window_assignment_matrix)
      window_basin_assign_sum <- window_basin_results$sum
      
      # Process HUC data for window production using custom function
      window_huc_result <- custom_process_huc_data(edges, basin, Huc, window_basin_assign_sum, HUC)
      
      # Extract necessary data from both results
      if (inherits(total_huc_result, "sf")) {
        total_huc_df <- st_drop_geometry(total_huc_result)
      } else {
        total_huc_df <- total_huc_result
      }
      
      if (inherits(window_huc_result, "sf")) {
        window_huc_df <- st_drop_geometry(window_huc_result)
      } else {
        window_huc_df <- window_huc_result
      }
      
      # Calculate total production for normalization
      total_production_all_hucs <- sum(total_huc_df$total_production, na.rm = TRUE)
      
      # Calculate the overall protection percentage for the entire run
      overall_protection_pct <- sum(window_huc_df$total_production, na.rm = TRUE) / 
        max(total_production_all_hucs, 0.0001) * 100  # Avoid division by zero
      
      # Join the datasets
      huc_col <- paste0("HUC", HUC)
      name_col <- "Name"
      
      # Make sure both dataframes have the necessary columns
      if (!huc_col %in% colnames(total_huc_df)) {
        message(paste("Warning: Missing", huc_col, "column in total_huc_df"))
        print(colnames(total_huc_df))
      }
      
      if (!name_col %in% colnames(total_huc_df)) {
        message(paste("Warning: Missing", name_col, "column in total_huc_df"))
        print(colnames(total_huc_df))
      }
      
      # Use a safer join approach
      combined_results <- total_huc_df %>%
        select(all_of(c(huc_col, name_col, "total_production"))) %>%
        rename(total_production_all = total_production) 
      
      # Add window production
      window_production <- window_huc_df %>%
        select(all_of(c(huc_col, "total_production"))) %>%
        rename(total_production_window = total_production)
      
      combined_results <- combined_results %>%
        left_join(window_production, by = huc_col)
      
      # Replace NA values with 0
      combined_results$total_production_window[is.na(combined_results$total_production_window)] <- 0
      
      # Calculate percentages
      combined_results <- combined_results %>%
        mutate(
          # Calculate total percentages (gray bars + green bars)
          huc_pct_of_total_run = (total_production_all / max(total_production_all_hucs, 0.0001)) * 100,
          
          # Calculate window percentages (green bars only)
          huc_pct_in_window = (total_production_window / max(total_production_all_hucs, 0.0001)) * 100,
          
          # Calculate the protection percentage (green as % of that HUC's total)
          protection_percentage = ifelse(total_production_all > 0,
                                         (total_production_window / total_production_all) * 100,
                                         0)
        ) %>%
        # Add year and metadata information
        mutate(
          year = year,
          watershed = watershed,
          closure_start_date = closure_start_date,
          closure_end_date = closure_end_date,
          closure_start_doy = closure_start_doy,
          closure_end_doy = closure_end_doy,
          overall_window_pct = overall_protection_pct
        ) %>%
        # Sort by HUC percentage of total run in descending order
        arrange(desc(huc_pct_of_total_run))
      
      # Add to overall results
      all_results <- bind_rows(all_results, combined_results)
      
    }, error = function(e) {
      message(paste("Error processing", watershed, "for year", year, ":", e$message))
    })
  }
  
  # Check if we have results
  if (nrow(all_results) == 0) {
    warning("No results collected - unable to save CSV")
    return(NULL)
  }
  
  # Create file name
  file_name <- paste0(
    "huc_closure_proportions_", 
    watershed, "_", 
    min(years), "_to_", max(years), "_", 
    closure_start_date, "_to_", closure_end_date, 
    ".csv"
  )
  
  output_path <- file.path(output_dir, file_name)
  
  # Save results to CSV
  write.csv(all_results, output_path, row.names = FALSE)
  
  message(paste("Saved HUC closure window proportions to:", output_path))
  
  # Create a summary file with just the key percentages
  summary_results <- all_results %>%
    select(year, Name, huc_pct_of_total_run, huc_pct_in_window, protection_percentage) %>%
    rename(
      "Year" = year,
      "HUC Name" = Name,
      "Total HUC Contribution %" = huc_pct_of_total_run,
      "Closure Window Contribution %" = huc_pct_in_window,
      "HUC Protection %" = protection_percentage
    )
  
  summary_file_name <- paste0(
    "huc_closure_summary_", 
    watershed, "_", 
    min(years), "_to_", max(years), "_", 
    closure_start_date, "_to_", closure_end_date, 
    ".csv"
  )
  
  summary_path <- file.path(output_dir, summary_file_name)
  write.csv(summary_results, summary_path, row.names = FALSE)
  
  message(paste("Saved simplified summary to:", summary_path))
  
  return(list(
    detailed = output_path,
    summary = summary_path
  ))
}

# Example usage
results <- export_huc_closure_proportions(
  years = c(2017, 2018, 2019, 2020, 2021),
  watershed = "Kusko",
  closure_start_date = "06-01", 
  closure_end_date = "06-11"
)

# To run for Yukon watershed, uncomment:
# results_yukon <- export_huc_closure_proportions(
#   years = c(2017, 2018, 2019, 2020, 2021),
#   watershed = "Yukon",
#   closure_start_date = "06-01", 
#   closure_end_date = "06-11"
# )