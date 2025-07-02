################################################################################
# MAIN.R - SIMPLIFIED SALMON RUN TIMING ANALYSIS (NO HUC MAPS)
################################################################################
# PURPOSE: Creates management river maps and exports CSV data:
#   1. DOY_Quartile: Shows proportion within each timing quartile
#   2. DOY_Total: Shows proportion of total annual run 
#   3. Front-end Closure: Shows protection during closure window
# NOTE: HUC map creation removed - only management river analysis
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

# Source required files (keeping existing structure)
source(here("code/utils/spatial_utils.R"))  # This now includes watershed ordering functions
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
# ANALYSIS 1: DOY QUARTILE ANALYSIS (MANAGEMENT RIVERS ONLY)
# PURPOSE: Shows what proportion of production occurs within each timing quartile
################################################################################

#' Run DOY Quartile Analysis with optional CSV export - Management Rivers Only
run_doy_analysis <- function(years, watersheds, export_csv = TRUE) {
  clean_outputs()
  
  # Storage for CSV export data - only management data
  all_mgmt_export_data <- NULL
  
  for (watershed in watersheds) {
    for (year in years) {
      params <- get_params(watershed)
      
      message(paste("Processing DOY quartiles (Management Rivers only):", year, watershed))
      
      # Run analysis and get return values for CSV export
      results <- DOY_Quartile_Analysis(
        year = year,
        watershed = watershed,
        sensitivity_threshold = params$sensitivity_threshold,
        min_error = params$min_error,
        min_stream_order = params$min_stream_order,
        return_values = export_csv  # Only get return values if we're exporting
      )
      
      # Process results for CSV export - only management data
      if (export_csv && !is.null(results)) {
        mgmt_export_data <- process_mgmt_results_for_export(results, year, watershed)
        
        # Store management data
        if (!is.null(mgmt_export_data)) {
          if (is.null(all_mgmt_export_data)) {
            all_mgmt_export_data <- mgmt_export_data
          } else {
            all_mgmt_export_data <- rbind(all_mgmt_export_data, mgmt_export_data)
          }
        }
      }
    }
  }
  
  # Export management CSV files if requested
  if (export_csv && !is.null(all_mgmt_export_data)) {
    export_mgmt_analysis_csv(all_mgmt_export_data)
  }
}

################################################################################
# ANALYSIS 2: DOY TOTAL ANALYSIS (MANAGEMENT RIVERS ONLY)
# PURPOSE: Shows what proportion of the TOTAL ANNUAL RUN each quartile represents
################################################################################

#' Run DOY Total Analysis with optional CSV export - Management Rivers Only
run_doy_total_analysis <- function(years, watersheds, export_csv = TRUE) {
  
  # Note: This analysis doesn't create HUC maps by default, only tributary maps
  # So we can keep it as-is but just skip any HUC-related processing
  
  for (watershed in watersheds) {
    for (year in years) {
      params <- get_params(watershed)
      
      message(paste("Processing DOY total (Management Rivers only):", year, watershed))
      
      # Run analysis - this will create tributary maps but not HUC maps
      results <- DOY_Total_Analysis(
        year = year,
        watershed = watershed,
        sensitivity_threshold = params$sensitivity_threshold,
        min_error = params$min_error,
        min_stream_order = params$min_stream_order,
        return_values = FALSE  # We don't need return values for this
      )
    }
  }
  
  message("DOY Total Analysis complete - tributary maps created")
}

################################################################################
# ANALYSIS 3: FRONT-END CLOSURE ANALYSIS
# PURPOSE: Shows what proportion of run is protected during closure window
################################################################################

#' Run Front-end Closure Analysis with optional CSV export
run_closure_analysis <- function(years, watershed = "Kusko", 
                                 closure_start = "06-01", closure_end = "06-11",
                                 export_csv = TRUE) {
  
  # Use the correct closure function from the right file
  source(here("Code/analysis/Front End Closure/FrontEndClosure.R"))
  
  # Run closure analysis and get results for CSV export
  closure_results <- analyze_closure_simple(
    years = years,
    watershed = watershed,
    closure_start = closure_start,
    closure_end = closure_end,
    create_maps = TRUE,
    clear_maps = TRUE
  )
  
  # Export CSV if requested and we have results
  if (export_csv && !is.null(closure_results)) {
    export_closure_analysis_csv(closure_results, closure_start, closure_end)
  }
  
  return(closure_results)
}

################################################################################
# CSV EXPORT FUNCTIONS FOR MANAGEMENT RIVER ANALYSIS
################################################################################

#' Process management river results for CSV export
process_mgmt_results_for_export <- function(results, year, watershed) {
  
  # Load natal data to calculate CPUE proportions
  natal_data <- load_natal_data(year, watershed)
  total_cpue <- sum(natal_data$dailyCPUEprop, na.rm = TRUE)
  
  # Divide into quartiles to get CPUE for each
  quartile_data <- divide_doy_quartiles(natal_data)
  quartile_cpue <- sapply(quartile_data$subsets, function(x) sum(x$dailyCPUEprop, na.rm = TRUE))
  quartile_cpue_prop <- quartile_cpue / total_cpue
  
  # Process each quartile result
  export_data <- NULL
  
  for (q in 1:length(results)) {
    if (is.null(results[[q]]) || is.null(results[[q]]$mgmt_result)) next
    
    mgmt_result <- results[[q]]$mgmt_result
    
    # Create export record for each management river
    quartile_data <- mgmt_result %>%
      select(mgmt_river, total_production, production_proportion, edge_count) %>%
      mutate(
        year = year,
        watershed = watershed,
        quartile = paste0("Q", q),
        quartile_num = q,
        within_quartile_prop = production_proportion,  # Proportion within this quartile
        cpue_prop_in_quartile = quartile_cpue_prop[q]  # CPUE proportion for this quartile
      ) %>%
      select(year, watershed, quartile, quartile_num, mgmt_river, 
             within_quartile_prop, cpue_prop_in_quartile, edge_count)
    
    export_data <- rbind(export_data, quartile_data)
  }
  
  return(export_data)
}

#' Calculate total run proportions for management river export data
calculate_mgmt_total_run_proportions <- function(export_data) {
  
  total_run_data <- NULL
  
  for (year in unique(export_data$year)) {
    for (watershed in unique(export_data$watershed)) {
      
      # Get parameters
      params <- get_params(watershed)
      
      # Load data
      spatial_data <- load_spatial_data(watershed, 8, params$min_stream_order)
      natal_data <- load_natal_data(year, watershed)
      quartile_data <- divide_doy_quartiles(natal_data)
      
      # Setup assignment parameters
      pid_iso <- spatial_data$edges$iso_pred
      pid_isose <- spatial_data$edges$isose_pred
      error <- calculate_error(pid_isose, params$min_error)
      priors <- setup_watershed_priors(spatial_data$edges, params$min_stream_order, watershed, natal_data)
      
      # Calculate total annual production
      all_data <- do.call(rbind, quartile_data$subsets)
      all_assignment_matrix <- perform_assignment(
        all_data, spatial_data$edges, watershed, priors, pid_iso, error, params$sensitivity_threshold
      )
      total_basin_assign_sum <- apply(all_assignment_matrix, 1, sum, na.rm = TRUE)
      grand_total_production <- sum(total_basin_assign_sum, na.rm = TRUE)
      
      # Process each quartile
      for (q in 1:4) {
        current_subset <- quartile_data$subsets[[q]]
        if (nrow(current_subset) == 0) next
        
        # Calculate quartile assignments
        assignment_matrix <- perform_assignment(
          current_subset, spatial_data$edges, watershed, priors, pid_iso, error, params$sensitivity_threshold
        )
        quartile_basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
        
        # Process management river data
        mgmt_result <- process_mgmt_river_data(spatial_data$edges, quartile_basin_assign_sum)
        
        if (!is.null(mgmt_result)) {
          # Calculate total run proportions
          total_run_props <- mgmt_result %>%
            mutate(
              year = year,
              watershed = watershed,
              quartile = paste0("Q", q),
              total_run_prop = total_production / grand_total_production
            ) %>%
            select(year, watershed, quartile, mgmt_river, total_run_prop)
          
          total_run_data <- rbind(total_run_data, total_run_props)
        }
      }
    }
  }
  
  return(total_run_data)
}

#' Export management river analysis data to CSV files
export_mgmt_analysis_csv <- function(export_data) {
  
  message("Calculating management river total run proportions...")
  total_run_data <- calculate_mgmt_total_run_proportions(export_data)
  
  # Merge export data with total run proportions
  final_data <- export_data %>%
    left_join(total_run_data, by = c("year", "watershed", "quartile", "mgmt_river"))
  
  # Create output directory
  csv_dir <- here("Analysis_Results/Management_River_Analysis")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  ################################################################################
  # EXPORT 1: TIDY FORMAT
  ################################################################################
  
  tidy_data <- final_data %>%
    select(
      year, 
      quartile,
      mgmt_river,
      total_run_prop,
      within_quartile_prop, 
      cpue_prop_in_quartile,
      edge_count
    ) %>%
    arrange(year, quartile, desc(total_run_prop))
  
  tidy_filepath <- file.path(csv_dir, "management_river_analysis_tidy.csv")
  write.csv(tidy_data, tidy_filepath, row.names = FALSE)
  message(paste("Exported management river tidy format to:", tidy_filepath))
  
  ################################################################################
  # EXPORT 2: TIMESERIES FORMAT - "Management River Timeseries.csv"
  ################################################################################
  
  # Create time period labels (sorted chronologically)
  final_data <- final_data %>%
    mutate(time_period = paste0(year, "_", quartile))
  
  # Get all unique time periods and sort them chronologically
  all_periods <- final_data %>%
    select(year, quartile_num, time_period) %>%
    distinct() %>%
    arrange(year, quartile_num) %>%
    pull(time_period)
  
  # Create timeseries for each management river with both metrics
  mgmt_timeseries <- NULL
  
  for (mgmt_name in unique(final_data$mgmt_river)) {
    mgmt_data <- final_data %>% filter(mgmt_river == mgmt_name)
    
    # Total run proportion row
    total_row <- mgmt_data %>%
      select(time_period, total_run_prop) %>%
      pivot_wider(names_from = time_period, values_from = total_run_prop, values_fill = 0) %>%
      mutate(
        Management_River = mgmt_name,
        Metric = "Total_Run_Proportion"
      ) %>%
      select(Management_River, Metric, all_of(all_periods))
    
    # Within quartile proportion row  
    within_row <- mgmt_data %>%
      select(time_period, within_quartile_prop) %>%
      pivot_wider(names_from = time_period, values_from = within_quartile_prop, values_fill = 0) %>%
      mutate(
        Management_River = mgmt_name,
        Metric = "Within_Quartile_Proportion"
      ) %>%
      select(Management_River, Metric, all_of(all_periods))
    
    # Combine both rows for this management river
    mgmt_timeseries <- bind_rows(mgmt_timeseries, total_row, within_row)
  }
  
  # Sort by management river name, then by metric
  mgmt_timeseries <- mgmt_timeseries %>%
    arrange(Management_River, Metric)
  
  timeseries_filepath <- file.path(csv_dir, "Management River Timeseries.csv")
  write.csv(mgmt_timeseries, timeseries_filepath, row.names = FALSE)
  message(paste("Exported management river timeseries to:", timeseries_filepath))
  
  # Print summary
  message("\n=== Management River CSV Export Summary ===")
  message(paste("Years processed:", paste(sort(unique(final_data$year)), collapse = ", ")))
  message(paste("Management rivers included:", length(unique(final_data$mgmt_river))))
  message(paste("Total records in tidy format:", nrow(tidy_data)))
  message(paste("Management timeseries rows:", nrow(mgmt_timeseries), "(", length(unique(final_data$mgmt_river)) * 2, "expected: 2 metrics per river)"))
}

################################################################################
# CSV EXPORT FUNCTIONS FOR CLOSURE ANALYSIS
################################################################################

#' Export Front-end Closure analysis data to CSV files
export_closure_analysis_csv <- function(closure_results, closure_start, closure_end) {
  
  # Create output directory
  csv_dir <- here("Analysis_Results/Closure_Analysis")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  ################################################################################
  # EXPORT 1: TIDY FORMAT
  ################################################################################
  
  # The analyze_closure_simple function returns columns:
  # Year, Management_Unit, Total_Contribution_Percent, Percent_Of_Unit_In_Window
  tidy_data <- closure_results %>%
    select(
      Year,
      Management_Unit,
      Total_Contribution_Percent,
      Percent_Of_Unit_In_Window
    ) %>%
    rename(
      year = Year,
      Name = Management_Unit,
      total_mgmt_contribution = Total_Contribution_Percent,
      closure_protection_percent = Percent_Of_Unit_In_Window
    ) %>%
    mutate(
      closure_window = paste0(closure_start, "_to_", closure_end),
      total_mgmt_contribution_prop = total_mgmt_contribution / 100,
      closure_protection_prop = closure_protection_percent / 100
    ) %>%
    arrange(year, desc(total_mgmt_contribution))
  
  tidy_filepath <- file.path(csv_dir, paste0("closure_analysis_tidy_", closure_start, "_to_", closure_end, ".csv"))
  write.csv(tidy_data, tidy_filepath, row.names = FALSE)
  message(paste("Exported Closure tidy format to:", tidy_filepath))
  
  ################################################################################
  # EXPORT 2: TIMESERIES FORMAT - "Closure Management Unit Timeseries.csv"
  ################################################################################
  
  # Create timeseries for each management unit
  mgmt_timeseries <- NULL
  
  for (mgmt_name in unique(tidy_data$Name)) {
    mgmt_data <- tidy_data %>% filter(Name == mgmt_name)
    
    # Total management unit contribution row
    total_row <- mgmt_data %>%
      select(year, total_mgmt_contribution_prop) %>%
      pivot_wider(names_from = year, values_from = total_mgmt_contribution_prop, values_fill = 0) %>%
      mutate(
        Management_Unit = mgmt_name,
        Metric = "Total_Management_Unit_Contribution"
      )
    
    # Closure protection row
    protection_row <- mgmt_data %>%
      select(year, closure_protection_prop) %>%
      pivot_wider(names_from = year, values_from = closure_protection_prop, values_fill = 0) %>%
      mutate(
        Management_Unit = mgmt_name,
        Metric = "Closure_Protection_Proportion"
      )
    
    # Combine both rows for this management unit
    mgmt_timeseries <- bind_rows(mgmt_timeseries, total_row, protection_row)
  }
  
  # Sort by management unit name, then by metric
  mgmt_timeseries <- mgmt_timeseries %>%
    arrange(Management_Unit, Metric)
  
  timeseries_filepath <- file.path(csv_dir, paste0("Closure Management Unit Timeseries_", closure_start, "_to_", closure_end, ".csv"))
  write.csv(mgmt_timeseries, timeseries_filepath, row.names = FALSE)
  message(paste("Exported Closure Management Unit timeseries to:", timeseries_filepath))
  
  message(paste("Closure CSV export complete -", nrow(tidy_data), "records,", length(unique(tidy_data$Name)), "management units"))
}

################################################################################
# MAIN EXECUTION SECTION
################################################################################

# Set analysis parameters
years <- c("2017","2018","2019","2020","2022")
watersheds <- c("Kusko")

message("=== Running Analysis Suite (Management Rivers Only) ===")

# 1. DOY Quartiles (Management Rivers only) + CSV export
message("Starting DOY Quartile Analysis (Management Rivers only) with CSV export...")
run_doy_analysis(years, watersheds, export_csv = TRUE)

# # 2. DOY Total (tributary maps only - no HUC maps)
# message("Starting DOY Total Analysis (Tributary maps only)...")
# run_doy_total_analysis(years, watersheds, export_csv = FALSE)

# 3. Front-end Closure + CSV export
message("Starting Front-end Closure Analysis with CSV export...")
run_closure_analysis(years, "Kusko", export_csv = TRUE)

message("=== Analysis Complete ===")
