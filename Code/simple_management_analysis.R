# simple_management_analysis.R
# Simple script that adds management unit analysis to your existing workflow
# This will work once you add the mgmt_river attribute to your original shapefile

library(sf)
library(dplyr)
library(here)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Clear environment
rm(list = ls())

cat("=== Simple Management Unit Analysis ===\n")
cat("This works with your original Kusko shapefile once mgmt_river attribute is added\n\n")

# Source the required files (using your existing functions)
source(here("code/utils/spatial_utils.R"))    # Reverted version
source(here("code/assignment.R"))              # Your existing assignment functions

#' Simple function to export management unit results alongside regular DOY analysis
#'
#' @param years Vector of years to process
#' @param output_dir Directory to save output CSV files
#' @return Path to management unit CSV file
export_management_unit_results <- function(years, output_dir = here("Analysis_Results")) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize data frame to store management results
  all_mgmt_data <- NULL
  
  # Parameters for Kusko
  watershed <- "Kusko"
  sensitivity_threshold <- 0.7
  min_error <- 0.0006
  min_stream_order <- 3
  
  # Process each year
  for (year in years) {
    cat("Processing year", year, "...\n")
    
    # Load spatial data
    spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    Huc <- spatial_data$Huc
    
    # Check if mgmt_river column exists
    if (!"mgmt_river" %in% colnames(edges)) {
      cat("  Warning: mgmt_river column not found for year", year, "- skipping\n")
      cat("  Please add mgmt_river attribute to your shapefile first\n")
      next
    }
    
    # Load natal origins data
    natal_data <- load_natal_data(year, watershed)
    
    # Divide into DOY quartiles
    quartile_data <- divide_doy_quartiles(natal_data)
    quartile_subsets <- quartile_data$subsets
    
    # Get assignment parameters
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, min_error)
    priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
    
    # Calculate total production across all quartiles
    all_data <- do.call(rbind, quartile_subsets)
    all_assignment_matrix <- perform_assignment(
      all_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    total_basin_assign_sum <- apply(all_assignment_matrix, 1, sum, na.rm = TRUE)
    grand_total_production <- sum(total_basin_assign_sum, na.rm = TRUE)
    
    cat("  Total production:", round(grand_total_production, 6), "\n")
    
    # Process each quartile
    for (q in 1:length(quartile_subsets)) {
      current_subset <- quartile_subsets[[q]]
      
      if (nrow(current_subset) == 0) next
      
      cat("  Processing Q", q, "with", nrow(current_subset), "data points\n")
      
      # Perform assignment for this quartile
      assignment_matrix <- perform_assignment(
        current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
      )
      quartile_basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
      
      # Process management unit data for this quartile
      mgmt_result <- process_management_data(edges, quartile_basin_assign_sum, grand_total_production)
      
      if (!is.null(mgmt_result) && nrow(mgmt_result) > 0) {
        # Add metadata
        mgmt_result$year <- year
        mgmt_result$quartile <- paste0("Q", q)
        mgmt_result$watershed <- watershed
        
        # Add to combined results
        all_mgmt_data <- rbind(all_mgmt_data, mgmt_result)
        
        cat("    Found", nrow(mgmt_result), "management units\n")
      } else {
        cat("    No management units found for Q", q, "\n")
      }
    }
  }
  
  # Save results
  if (!is.null(all_mgmt_data) && nrow(all_mgmt_data) > 0) {
    output_file <- file.path(output_dir, "management_unit_analysis.csv")
    write.csv(all_mgmt_data, output_file, row.names = FALSE)
    cat("\nSaved management unit results to:", output_file, "\n")
    
    # Show summary
    cat("\nSummary of results:\n")
    summary_table <- all_mgmt_data %>%
      group_by(mgmt_river) %>%
      summarise(
        total_contribution_percent = sum(percent_of_total_run),
        avg_per_quartile = mean(percent_of_total_run),
        max_contribution = max(percent_of_total_run),
        years_present = n_distinct(year),
        .groups = "drop"
      ) %>%
      arrange(desc(total_contribution_percent))
    
    print(summary_table)
    
    return(output_file)
  } else {
    cat("\nNo management unit data collected.\n")
    cat("Make sure your shapefile has the mgmt_river attribute.\n")
    return(NULL)
  }
}

#' Test if management unit analysis will work
test_management_setup <- function() {
  cat("Testing management unit setup...\n")
  
  # Try to load one year of data
  tryCatch({
    spatial_data <- load_spatial_data("Kusko", 8, 3)
    edges <- spatial_data$edges
    
    cat("Successfully loaded", nrow(edges), "edges\n")
    
    if ("mgmt_river" %in% colnames(edges)) {
      managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
      cat("Found", nrow(managed_edges), "edges with management unit assignments\n")
      
      if (nrow(managed_edges) > 0) {
        unique_mgmt <- unique(managed_edges$mgmt_river)
        cat("Management rivers:", paste(unique_mgmt, collapse = ", "), "\n")
        cat("Setup looks good! Ready to run analysis.\n")
        return(TRUE)
      } else {
        cat("No edges have management unit assignments\n")
        return(FALSE)
      }
    } else {
      cat("mgmt_river column not found\n")
      cat("You need to add this attribute to your shapefile first\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat("Error testing setup:", e$message, "\n")
    return(FALSE)
  })
}

# Run the test
cat("Step 1: Testing if your shapefile is ready...\n")
ready <- test_management_setup()

if (ready) {
  cat("\nStep 2: Running management unit analysis...\n")
  years_to_analyze <- c("2017", "2019", "2020", "2021")
  
  result_file <- export_management_unit_results(years_to_analyze)
  
  if (!is.null(result_file)) {
    cat("\n=== Analysis Complete! ===\n")
    cat("Management unit data saved to:", result_file, "\n")
    cat("The 'percent_of_total_run' column shows what % of the entire watershed run each management unit represents\n")
  }
} else {
  cat("\n=== Setup Not Ready ===\n")
  cat("Please add the mgmt_river attribute to your original Kusko shapefile first.\n")
  cat("Once you do that, run this script again.\n")
}