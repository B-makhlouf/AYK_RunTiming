################################################################################
# 01_TRIBUTARY_MAPS.R - DOY QUARTILE TRIBUTARY MAPPING ANALYSIS
################################################################################
# Creates tributary maps showing spatial distribution of production by timing quartile
# Run 00_setup.R and 05_visualization_functions.R first
################################################################################

cat("=== TRIBUTARY MAPPING ANALYSIS ===\n")

# Check if setup is loaded
if (!exists("CONFIG")) {
  stop("Please run 00_setup.R first to load parameters and functions")
}

if (!exists("create_tributary_map")) {
  stop("Please run 05_visualization_functions.R first to load plotting functions")
}

################################################################################
# MAIN TRIBUTARY ANALYSIS FUNCTION
################################################################################

#' Run tributary mapping analysis for DOY quartiles
#' Creates maps showing spatial distribution of production by timing quartile
#' Produces identical outputs to original DOY_analysis.R
run_tributary_analysis <- function(years = CONFIG$years, 
                                   watersheds = CONFIG$watersheds,
                                   export_csv = TRUE) {
  
  cat("Starting tributary mapping analysis...\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n")
  
  # Create output directories
  create_output_dirs()
  tributary_dir <- file.path(PATHS$maps_dir, "DOY_Quartile", "Tribs")
  mgmt_dir <- file.path(PATHS$maps_dir, "DOY_Quartile", "Management")
  dir.create(tributary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mgmt_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Storage for CSV export data
  all_mgmt_export_data <- NULL
  
  for (watershed in watersheds) {
    params <- WATERSHED_PARAMS[[watershed]]
    
    for (year in years) {
      cat(glue("Processing {watershed} {year}...\n"))
      
      # Load data
      spatial_data <- load_spatial_data(watershed)
      natal_data <- load_natal_data(year, watershed)
      
      cat(glue("  Loaded {nrow(natal_data)} fish observations\n"))
      
      # Setup for assignment
      pid_iso <- spatial_data$edges$iso_pred
      pid_isose <- spatial_data$edges$isose_pred
      error <- calculate_error(pid_isose, params$min_error)
      priors <- setup_priors(spatial_data$edges, watershed, natal_data)
      
      # Divide into quartiles
      quartile_data <- divide_doy_quartiles(natal_data)
      
      # Storage for year's results
      year_results <- list()
      
      # Process each quartile
      for (q in 1:4) {
        current_subset <- quartile_data$subsets[[q]]
        
        if (nrow(current_subset) == 0) {
          cat(glue("  Skipping Q{q} - no data\n"))
          next
        }
        
        cat(glue("  Processing Q{q} with {nrow(current_subset)} fish\n"))
        
        #----------------------------------------------------------------------
        # PERFORM ASSIGNMENT
        #----------------------------------------------------------------------
        
        assignment_matrix <- perform_assignment(
          current_subset, spatial_data$edges, watershed, priors, 
          pid_iso, error, params$sensitivity_threshold
        )
        
        # Calculate basin assignments
        basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
        basin_assign_rescale <- basin_assign_sum / sum(basin_assign_sum, na.rm = TRUE)
        basin_assign_norm <- basin_assign_rescale / max(basin_assign_rescale, na.rm = TRUE)
        
        #----------------------------------------------------------------------
        # CREATE TRIBUTARY MAP WITH HISTOGRAM
        #----------------------------------------------------------------------
        
        # Create DOY histogram
        gg_hist <- create_doy_histogram(natal_data, current_subset, quartile_data$labels[q])
        
        # Create tributary map (matching original style exactly) - NOW WITH HISTOGRAM
        trib_filename <- glue("{watershed}_{year}_DOY_Q{q}.png")
        trib_path <- file.path(tributary_dir, trib_filename)
        
        # [trimmed: not a paper figure] per-year/quartile tributary production
        # maps are not exported. The management-unit aggregation below still runs
        # and writes management_river_analysis_tidy.csv (needed by step 06).
        # create_tributary_map(
        #   basin = spatial_data$basin,
        #   edges = spatial_data$edges,
        #   basin_assign_norm = basin_assign_norm,
        #   year = year,
        #   watershed = watershed,
        #   subset_label = quartile_data$labels[q],
        #   output_path = trib_path,
        #   StreamOrderPrior = priors$StreamOrderPrior,
        #   pid_prior = priors$pid_prior,
        #   gg_hist = gg_hist
        # )
        
        #----------------------------------------------------------------------
        # PROCESS MANAGEMENT RIVER DATA
        #----------------------------------------------------------------------
        
        mgmt_result <- process_mgmt_assignments(spatial_data$edges, basin_assign_rescale)
        
        if (!is.null(mgmt_result)) {
          # Create management river map
          mgmt_filename <- glue("{watershed}_{year}_DOY_Q{q}_Management.png")
          mgmt_path <- file.path(mgmt_dir, mgmt_filename)
          
          # [trimmed: not a paper figure] per-year/quartile management maps
          # are not exported. Aggregation/CSV export below is unchanged.
          # create_mgmt_river_map(
          #   mgmt_result, spatial_data$edges, spatial_data$basin,
          #   year, watershed, quartile_data$labels[q], mgmt_path, gg_hist
          # )
          
          # Store results for CSV export
          year_results[[q]] <- list(
            subset = current_subset,
            label = quartile_data$labels[q],
            basin_assign_rescale = basin_assign_rescale,
            basin_assign_norm = basin_assign_norm,
            mgmt_result = mgmt_result
          )
        }
      }
      
      #------------------------------------------------------------------------
      # PROCESS RESULTS FOR CSV EXPORT
      #------------------------------------------------------------------------
      
      if (export_csv && length(year_results) > 0) {
        mgmt_export_data <- process_mgmt_results_for_export(year_results, year, watershed)
        
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
  
  #--------------------------------------------------------------------------
  # EXPORT CSV FILES
  #--------------------------------------------------------------------------
  
  if (export_csv && !is.null(all_mgmt_export_data)) {
    export_mgmt_analysis_csv(all_mgmt_export_data)
  }
  
  cat("✓ Tributary mapping analysis complete\n")
  cat("  Maps saved to:\n")
  cat(glue("    - Tributary: {tributary_dir}\n"))
  cat(glue("    - Management: {mgmt_dir}\n"))
  if (export_csv) {
    cat("  CSV data exported to Analysis_Results/Management_River_Analysis/\n")
  }
}

################################################################################
# CSV EXPORT FUNCTIONS (matching original exactly)
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
      params <- WATERSHED_PARAMS[[watershed]]
      
      # Load data
      spatial_data <- load_spatial_data(watershed)
      natal_data <- load_natal_data(year, watershed)
      quartile_data <- divide_doy_quartiles(natal_data)
      
      # Setup assignment parameters
      pid_iso <- spatial_data$edges$iso_pred
      pid_isose <- spatial_data$edges$isose_pred
      error <- calculate_error(pid_isose, params$min_error)
      priors <- setup_priors(spatial_data$edges, watershed, natal_data)
      
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
        mgmt_result <- process_mgmt_assignments(spatial_data$edges, quartile_basin_assign_sum)
        
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

#' Export management river analysis data to CSV files (matching original exactly)
export_mgmt_analysis_csv <- function(export_data) {
  
  cat("Calculating management river total run proportions...\n")
  total_run_data <- calculate_mgmt_total_run_proportions(export_data)
  
  # Merge export data with total run proportions
  final_data <- export_data %>%
    left_join(total_run_data, by = c("year", "watershed", "quartile", "mgmt_river"))
  
  # Create output directory
  csv_dir <- file.path(PATHS$output_dir, "Management_River_Analysis")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  #--------------------------------------------------------------------------
  # EXPORT 1: TIDY FORMAT
  #--------------------------------------------------------------------------
  
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
  write_csv(tidy_data, tidy_filepath)
  cat(glue("Exported management river tidy format to: {tidy_filepath}\n"))
  
  #--------------------------------------------------------------------------
  # EXPORT 2: TIMESERIES FORMAT
  #--------------------------------------------------------------------------
  
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
  write_csv(mgmt_timeseries, timeseries_filepath)
  cat(glue("Exported management river timeseries to: {timeseries_filepath}\n"))
  
  # Print summary
  cat("\n=== Management River CSV Export Summary ===\n")
  cat(glue("Years processed: {paste(sort(unique(final_data$year)), collapse = ', ')}\n"))
  cat(glue("Management rivers included: {length(unique(final_data$mgmt_river))}\n"))
  cat(glue("Total records in tidy format: {nrow(tidy_data)}\n"))
  cat(glue("Management timeseries rows: {nrow(mgmt_timeseries)} ({length(unique(final_data$mgmt_river)) * 2} expected: 2 metrics per river)\n"))
}

################################################################################
# EXECUTION SECTION
################################################################################

# Run the analysis if this script is executed directly
if (interactive() || !exists(".script_executed")) {
  cat("Running tributary mapping analysis...\n")
  cat("Make sure you have run:\n")
  cat("  source('00_setup.R')\n")
  cat("  source('05_visualization_functions.R')\n")
  cat("before running this script.\n\n")
  
  # Uncomment the line below to run the analysis
  # run_tributary_analysis()
  
  .script_executed <- TRUE
}

cat("✓ Tributary mapping functions loaded\n")
cat("Run: run_tributary_analysis() to execute the analysis\n")