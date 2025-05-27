################################################################################
# MAIN.R - SIMPLIFIED SALMON RUN TIMING ANALYSIS
################################################################################
# PURPOSE: Creates three types of maps and exports CSV data:
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
library(tidyr)

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

#' Run DOY Quartile Analysis with optional CSV export
run_doy_analysis <- function(years, watersheds, export_csv = TRUE) {
  clean_outputs()
  
  # Storage for CSV export data
  all_export_data <- NULL
  all_mgmt_export_data <- NULL  # Add management data storage
  
  for (watershed in watersheds) {
    for (year in years) {
      params <- get_params(watershed)
      
      message(paste("Processing DOY quartiles:", year, watershed))
      
      # Run analysis and get return values for CSV export
      results <- DOY_Quartile_Analysis(
        year = year,
        watershed = watershed,
        sensitivity_threshold = params$sensitivity_threshold,
        min_error = params$min_error,
        min_stream_order = params$min_stream_order,
        HUC = 8,
        return_values = export_csv  # Only get return values if we're exporting
      )
      
      # Process results for CSV export
      if (export_csv && !is.null(results)) {
        export_data <- process_results_for_export(results, year, watershed)
        mgmt_export_data <- process_mgmt_results_for_export(results, year, watershed)
        all_export_data <- rbind(all_export_data, export_data)
        
        # Store management data separately
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
  
  # Export CSV files if requested
  if (export_csv && !is.null(all_export_data)) {
    export_doy_analysis_csv(all_export_data)
  }
  
  # Export management CSV files if requested
  if (export_csv && !is.null(all_mgmt_export_data)) {
    export_mgmt_analysis_csv(all_mgmt_export_data)
  }
}

################################################################################
# ANALYSIS 2: DOY TOTAL ANALYSIS  
# PURPOSE: Shows what proportion of the TOTAL ANNUAL RUN each quartile represents
################################################################################

#' Run DOY Total Analysis with optional CSV export
run_doy_total_analysis <- function(years, watersheds, export_csv = TRUE) {
  
  # Storage for CSV export data
  all_total_export_data <- NULL
  
  for (watershed in watersheds) {
    for (year in years) {
      params <- get_params(watershed)
      
      message(paste("Processing DOY total:", year, watershed))
      
      # Run analysis and get return values for CSV export
      results <- DOY_Total_Analysis(
        year = year,
        watershed = watershed,
        sensitivity_threshold = params$sensitivity_threshold,
        min_error = params$min_error,
        min_stream_order = params$min_stream_order,
        HUC = 8,
        return_values = export_csv  # Only get return values if we're exporting
      )
      
      # Process results for CSV export
      if (export_csv && !is.null(results)) {
        export_data <- process_total_results_for_export(results, year, watershed)
        all_total_export_data <- rbind(all_total_export_data, export_data)
      }
    }
  }
  
  # Export CSV files if requested
  if (export_csv && !is.null(all_total_export_data)) {
    export_doy_total_analysis_csv(all_total_export_data)
  }
}

################################################################################
# ANALYSIS 3: FRONT-END CLOSURE ANALYSIS
# PURPOSE: Shows what proportion of run is protected during closure window
################################################################################

#' Run Front-end Closure Analysis with optional CSV export
run_closure_analysis <- function(years, watershed = "Kusko", 
                                 closure_start = "06-01", closure_end = "06-11",
                                 export_csv = TRUE) {
  
  # Use existing closure function from FrontEndClosure.R
  source(here("code/FrontEndClosure.R"))
  
  # Run closure analysis and get results for CSV export
  closure_results <- analyze_closure_window(
    years = years,
    watershed = watershed,
    closure_start = closure_start,
    closure_end = closure_end,
    create_maps = TRUE
  )
  
  # Export CSV if requested and we have results
  if (export_csv && !is.null(closure_results)) {
    export_closure_analysis_csv(closure_results, closure_start, closure_end)
  }
}

################################################################################
# CSV EXPORT FUNCTIONS
################################################################################

#' Process DOY analysis results for CSV export
process_results_for_export <- function(results, year, watershed) {
  
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
    if (is.null(results[[q]])) next
    
    huc_result <- results[[q]]$huc_result
    
    # Extract HUC data
    if (inherits(huc_result, "sf")) {
      huc_df <- st_drop_geometry(huc_result)
    } else {
      huc_df <- huc_result
    }
    
    # Create export record for each HUC
    quartile_data <- huc_df %>%
      select(Name, total_production, production_proportion) %>%
      mutate(
        year = year,
        watershed = watershed,
        quartile = paste0("Q", q),
        quartile_num = q,
        within_quartile_prop = production_proportion,  # Proportion within this quartile
        cpue_prop_in_quartile = quartile_cpue_prop[q]  # CPUE proportion for this quartile
      ) %>%
      select(year, watershed, quartile, quartile_num, Name, 
             within_quartile_prop, cpue_prop_in_quartile)
    
    export_data <- rbind(export_data, quartile_data)
  }
  
  return(export_data)
}

#' Calculate total run proportions for export data
calculate_total_run_proportions <- function(export_data) {
  
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
        
        # Process HUC data
        huc_result <- process_huc_data(spatial_data$edges, spatial_data$basin, 
                                       spatial_data$Huc, quartile_basin_assign_sum, 8)
        
        # Calculate total run proportions
        huc_df <- if (inherits(huc_result, "sf")) st_drop_geometry(huc_result) else huc_result
        
        total_run_props <- huc_df %>%
          mutate(
            year = year,
            watershed = watershed,
            quartile = paste0("Q", q),
            total_run_prop = total_production / grand_total_production
          ) %>%
          select(year, watershed, quartile, Name, total_run_prop)
        
        total_run_data <- rbind(total_run_data, total_run_props)
      }
    }
  }
  
  return(total_run_data)
}

#' Export DOY analysis data to CSV files
export_doy_analysis_csv <- function(export_data) {
  
  message("Calculating total run proportions...")
  total_run_data <- calculate_total_run_proportions(export_data)
  
  # Merge export data with total run proportions
  final_data <- export_data %>%
    left_join(total_run_data, by = c("year", "watershed", "quartile", "Name"))
  
  # Create output directory
  csv_dir <- here("Analysis_Results/DOY_Analysis")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  ################################################################################
  # EXPORT 1: TIDY FORMAT
  ################################################################################
  
  tidy_data <- final_data %>%
    select(
      year, 
      quartile,
      Name,
      total_run_prop,
      within_quartile_prop, 
      cpue_prop_in_quartile
    ) %>%
    arrange(year, quartile, desc(total_run_prop))
  
  tidy_filepath <- file.path(csv_dir, "doy_quartile_analysis_tidy.csv")
  write.csv(tidy_data, tidy_filepath, row.names = FALSE)
  message(paste("Exported tidy format to:", tidy_filepath))
  
  ################################################################################
  # EXPORT 2: TIMESERIES FORMAT - "HUC Timeseries.csv"
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
  
  # Create timeseries for each HUC with both metrics
  huc_timeseries <- NULL
  
  for (huc_name in unique(final_data$Name)) {
    huc_data <- final_data %>% filter(Name == huc_name)
    
    # Total run proportion row
    total_row <- huc_data %>%
      select(time_period, total_run_prop) %>%
      pivot_wider(names_from = time_period, values_from = total_run_prop, values_fill = 0) %>%
      mutate(
        HUC = huc_name,
        Metric = "Total_Run_Proportion"
      ) %>%
      select(HUC, Metric, all_of(all_periods))
    
    # Within quartile proportion row  
    within_row <- huc_data %>%
      select(time_period, within_quartile_prop) %>%
      pivot_wider(names_from = time_period, values_from = within_quartile_prop, values_fill = 0) %>%
      mutate(
        HUC = huc_name,
        Metric = "Within_Quartile_Proportion"
      ) %>%
      select(HUC, Metric, all_of(all_periods))
    
    # Combine both rows for this HUC
    huc_timeseries <- bind_rows(huc_timeseries, total_row, within_row)
  }
  
  # Sort by HUC name, then by metric (Total_Run first, then Within_Quartile)
  huc_timeseries <- huc_timeseries %>%
    arrange(HUC, Metric)
  
  timeseries_filepath <- file.path(csv_dir, "HUC Timeseries.csv")
  write.csv(huc_timeseries, timeseries_filepath, row.names = FALSE)
  message(paste("Exported HUC timeseries to:", timeseries_filepath))
  
  # Print summary
  message("\n=== CSV Export Summary ===")
  message(paste("Years processed:", paste(sort(unique(final_data$year)), collapse = ", ")))
  message(paste("HUCs included:", length(unique(final_data$Name))))
  message(paste("Total records in tidy format:", nrow(tidy_data)))
  message(paste("HUC timeseries rows:", nrow(huc_timeseries), "(", length(unique(final_data$Name)) * 2, "expected: 2 metrics per HUC)"))
}

################################################################################
# CSV EXPORT FUNCTIONS FOR DOY TOTAL ANALYSIS
################################################################################

#' Process DOY Total analysis results for CSV export
process_total_results_for_export <- function(results, year, watershed) {
  
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
    if (is.null(results[[q]])) next
    
    huc_result <- results[[q]]$huc_result
    
    # Extract HUC data
    if (inherits(huc_result, "sf")) {
      huc_df <- st_drop_geometry(huc_result)
    } else {
      huc_df <- huc_result
    }
    
    # Create export record for each HUC
    quartile_data <- huc_df %>%
      select(Name, total_production, percent_of_total_run) %>%
      mutate(
        year = year,
        watershed = watershed,
        quartile = paste0("Q", q),
        quartile_num = q,
        total_run_prop = percent_of_total_run / 100,  # Convert from percentage
        cpue_prop_in_quartile = quartile_cpue_prop[q]  # CPUE proportion for this quartile
      ) %>%
      select(year, watershed, quartile, quartile_num, Name, 
             total_run_prop, cpue_prop_in_quartile)
    
    export_data <- rbind(export_data, quartile_data)
  }
  
  return(export_data)
}

#' Export DOY Total analysis data to CSV files
export_doy_total_analysis_csv <- function(export_data) {
  
  # Create output directory
  csv_dir <- here("Analysis_Results/DOY_Total_Analysis")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  ################################################################################
  # EXPORT 1: TIDY FORMAT
  ################################################################################
  
  tidy_data <- export_data %>%
    select(
      year, 
      quartile,
      Name,
      total_run_prop,
      cpue_prop_in_quartile
    ) %>%
    arrange(year, quartile, desc(total_run_prop))
  
  tidy_filepath <- file.path(csv_dir, "doy_total_analysis_tidy.csv")
  write.csv(tidy_data, tidy_filepath, row.names = FALSE)
  message(paste("Exported DOY Total tidy format to:", tidy_filepath))
  
  ################################################################################
  # EXPORT 2: TIMESERIES FORMAT - "DOY Total HUC Timeseries.csv"
  ################################################################################
  
  # Create time period labels (sorted chronologically)
  export_data <- export_data %>%
    mutate(time_period = paste0(year, "_", quartile))
  
  # Get all unique time periods and sort them chronologically
  all_periods <- export_data %>%
    select(year, quartile_num, time_period) %>%
    distinct() %>%
    arrange(year, quartile_num) %>%
    pull(time_period)
  
  # Create timeseries for each HUC
  huc_timeseries <- NULL
  
  for (huc_name in unique(export_data$Name)) {
    huc_data <- export_data %>% filter(Name == huc_name)
    
    # Total run proportion row
    total_row <- huc_data %>%
      select(time_period, total_run_prop) %>%
      pivot_wider(names_from = time_period, values_from = total_run_prop, values_fill = 0) %>%
      mutate(
        HUC = huc_name,
        Metric = "Total_Run_Proportion"
      ) %>%
      select(HUC, Metric, all_of(all_periods))
    
    huc_timeseries <- bind_rows(huc_timeseries, total_row)
  }
  
  # Sort by HUC name
  huc_timeseries <- huc_timeseries %>%
    arrange(HUC)
  
  timeseries_filepath <- file.path(csv_dir, "DOY Total HUC Timeseries.csv")
  write.csv(huc_timeseries, timeseries_filepath, row.names = FALSE)
  message(paste("Exported DOY Total HUC timeseries to:", timeseries_filepath))
  
  message(paste("DOY Total CSV export complete -", nrow(tidy_data), "records,", length(unique(export_data$Name)), "HUCs"))
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
  
  tidy_data <- closure_results %>%
    select(
      Year,
      HUC_Name,
      Total_HUC_Contribution_Percent,
      Proportion_in_Closure_Window_Percent
    ) %>%
    rename(
      year = Year,
      Name = HUC_Name,
      total_huc_contribution = Total_HUC_Contribution_Percent,
      closure_protection_percent = Proportion_in_Closure_Window_Percent
    ) %>%
    mutate(
      closure_window = paste0(closure_start, "_to_", closure_end),
      total_huc_contribution_prop = total_huc_contribution / 100,
      closure_protection_prop = closure_protection_percent / 100
    ) %>%
    arrange(year, desc(total_huc_contribution))
  
  tidy_filepath <- file.path(csv_dir, paste0("closure_analysis_tidy_", closure_start, "_to_", closure_end, ".csv"))
  write.csv(tidy_data, tidy_filepath, row.names = FALSE)
  message(paste("Exported Closure tidy format to:", tidy_filepath))
  
  ################################################################################
  # EXPORT 2: TIMESERIES FORMAT - "Closure HUC Timeseries.csv"
  ################################################################################
  
  # Create timeseries for each HUC
  huc_timeseries <- NULL
  
  for (huc_name in unique(tidy_data$Name)) {
    huc_data <- tidy_data %>% filter(Name == huc_name)
    
    # Total HUC contribution row
    total_row <- huc_data %>%
      select(year, total_huc_contribution_prop) %>%
      pivot_wider(names_from = year, values_from = total_huc_contribution_prop, values_fill = 0) %>%
      mutate(
        HUC = huc_name,
        Metric = "Total_HUC_Contribution"
      )
    
    # Closure protection row
    protection_row <- huc_data %>%
      select(year, closure_protection_prop) %>%
      pivot_wider(names_from = year, values_from = closure_protection_prop, values_fill = 0) %>%
      mutate(
        HUC = huc_name,
        Metric = "Closure_Protection_Proportion"
      )
    
    # Combine both rows for this HUC
    huc_timeseries <- bind_rows(huc_timeseries, total_row, protection_row)
  }
  
  # Sort by HUC name, then by metric
  huc_timeseries <- huc_timeseries %>%
    arrange(HUC, Metric)
  
  timeseries_filepath <- file.path(csv_dir, paste0("Closure HUC Timeseries_", closure_start, "_to_", closure_end, ".csv"))
  write.csv(huc_timeseries, timeseries_filepath, row.names = FALSE)
  message(paste("Exported Closure HUC timeseries to:", timeseries_filepath))
  
  message(paste("Closure CSV export complete -", nrow(tidy_data), "records,", length(unique(tidy_data$Name)), "HUCs"))
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
# MAIN EXECUTION SECTION
################################################################################

# Set analysis parameters
years <- c("2017", "2019", "2020", "2021")
watersheds <- c("Kusko")

message("=== Running Simplified Analysis Suite ===")

# 1. DOY Quartiles (proportion within each quartile) + CSV export
message("Starting DOY Quartile Analysis with CSV export...")
run_doy_analysis(years, watersheds, export_csv = TRUE)

# 2. DOY Total (proportion of total run) + CSV export
message("Starting DOY Total Analysis with CSV export...")
run_doy_total_analysis(years, watersheds, export_csv = TRUE)

# 3. Front-end Closure + CSV export
message("Starting Front-end Closur Analysis with CSV export...")
run_closure_analysis(years, "Kusko", export_csv = TRUE)

message("=== Analysis Complete ===")