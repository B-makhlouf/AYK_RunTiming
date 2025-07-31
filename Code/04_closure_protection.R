################################################################################
# 04_CLOSURE_PROTECTION.R - FRONT-END CLOSURE PROTECTION ANALYSIS
################################################################################
# Analyzes what percentage of production is protected during closure periods
# Creates protection effectiveness boxplots and summary statistics
# Run 00_setup.R and 05_visualization_functions.R first
################################################################################

cat("=== FRONT-END CLOSURE PROTECTION ANALYSIS ===\n")

# Check if setup is loaded
if (!exists("CONFIG")) {
  stop("Please run 00_setup.R first to load parameters and functions")
}

if (!exists("create_closure_boxplot")) {
  stop("Please run 05_visualization_functions.R first to load plotting functions")
}

################################################################################
# MAIN CLOSURE PROTECTION ANALYSIS FUNCTION
################################################################################

#' Run front-end closure protection analysis  
#' Analyzes what percentage of production is protected during closure periods
#' Produces identical outputs to original closure protection analysis
run_closure_analysis <- function(years = CONFIG$years,
                                 watersheds = CONFIG$watersheds,
                                 closure_start = CONFIG$closure_dates$start,
                                 closure_end = CONFIG$closure_dates$end,
                                 export_csv = TRUE) {
  
  cat("Starting front-end closure protection analysis...\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n")
  cat("Closure window:", closure_start, "to", closure_end, "\n")
  
  # Create output directory
  output_dir <- file.path(PATHS$output_dir, "Closure_Protection")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Convert closure dates to DOY
  start_doy <- as.numeric(format(as.Date(paste0("2024-", closure_start)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2024-", closure_end)), "%j"))
  
  cat(glue("Analyzing closure window: DOY {start_doy} to {end_doy}\n"))
  
  all_results <- list()
  
  for (watershed in watersheds) {
    params <- WATERSHED_PARAMS[[watershed]]
    
    for (year in years) {
      cat(glue("Processing closure analysis for {watershed} {year}...\n"))
      
      # Load data
      spatial_data <- load_spatial_data(watershed)
      natal_data <- load_natal_data(year, watershed)
      
      cat(glue("  Loaded {nrow(natal_data)} fish observations\n"))
      
      # Check for management units
      if (!"mgmt_river" %in% colnames(spatial_data$edges)) {
        warning("Management river data not found. Skipping year ", year)
        next
      }
      
      # Setup assignment parameters
      pid_iso <- spatial_data$edges$iso_pred
      pid_isose <- spatial_data$edges$isose_pred
      error <- calculate_error(pid_isose, params$min_error)
      priors <- setup_priors(spatial_data$edges, watershed, natal_data)
      
      #------------------------------------------------------------------------
      # CALCULATE TOTAL ANNUAL PRODUCTION
      #------------------------------------------------------------------------
      
      cat("  Calculating total annual production...\n")
      total_assignment <- perform_assignment(
        natal_data, spatial_data$edges, watershed, priors,
        pid_iso, error, params$sensitivity_threshold
      )
      total_basin_assign <- apply(total_assignment, 1, sum, na.rm = TRUE)
      total_mgmt_result <- process_mgmt_assignments(spatial_data$edges, total_basin_assign)
      
      #------------------------------------------------------------------------
      # CALCULATE CLOSURE WINDOW PRODUCTION
      #------------------------------------------------------------------------
      
      closure_data <- natal_data %>% filter(DOY >= start_doy, DOY <= end_doy)
      
      cat(glue("  Found {nrow(closure_data)} fish in closure window\n"))
      
      if (nrow(closure_data) > 0 && !is.null(total_mgmt_result)) {
        closure_assignment <- perform_assignment(
          closure_data, spatial_data$edges, watershed, priors,
          pid_iso, error, params$sensitivity_threshold
        )
        closure_basin_assign <- apply(closure_assignment, 1, sum, na.rm = TRUE)
        closure_mgmt_result <- process_mgmt_assignments(spatial_data$edges, closure_basin_assign)
        
        #----------------------------------------------------------------------
        # CALCULATE PROTECTION PERCENTAGES
        #----------------------------------------------------------------------
        
        if (!is.null(closure_mgmt_result)) {
          protection_results <- total_mgmt_result %>%
            left_join(closure_mgmt_result %>% 
                        select(mgmt_river, closure_production = total_production),
                      by = "mgmt_river") %>%
            mutate(
              year = year,
              watershed = watershed,
              closure_production = ifelse(is.na(closure_production), 0, closure_production),
              protection_percentage = (closure_production / total_production) * 100,
              closure_window = glue("{closure_start}_to_{closure_end}"),
              total_contribution_percent = (total_production / sum(total_production)) * 100
            )
          
          all_results[[glue("{watershed}_{year}")]] <- protection_results
          
          cat(glue("  Processed {nrow(protection_results)} management units\n"))
        }
      }
    }
  }
  
  #--------------------------------------------------------------------------
  # COMBINE AND PROCESS RESULTS
  #--------------------------------------------------------------------------
  
  if (length(all_results) == 0) {
    stop("No valid results processed for any year")
  }
  
  # Combine results
  combined_results <- bind_rows(all_results)
  
  cat(glue("Combined results: {nrow(combined_results)} records\n"))
  cat(glue("Management units: {length(unique(combined_results$mgmt_river))}\n"))
  
  # Apply watershed ordering for consistency (reverse for boxplots)
  combined_results <- apply_watershed_order(combined_results, reverse_for_plots = TRUE)
  
  # Print summary statistics
  overall_protection <- combined_results %>%
    group_by(year) %>%
    summarise(
      total_annual_production = sum(total_production),
      total_closure_production = sum(closure_production),
      overall_protection_pct = (total_closure_production / total_annual_production) * 100,
      .groups = "drop"
    )
  
  cat("\nOverall protection by year:\n")
  print(overall_protection)
  
  #--------------------------------------------------------------------------
  # CREATE VISUALIZATIONS
  #--------------------------------------------------------------------------
  
  create_closure_plots(combined_results, output_dir)  # Creates BOTH plots
  
  
  #--------------------------------------------------------------------------
  # CREATE SUMMARY STATISTICS
  #--------------------------------------------------------------------------
  
  # Calculate summary statistics by management unit
  summary_stats <- combined_results %>%
    group_by(mgmt_river) %>%
    summarise(
      mean_protection_pct = round(mean(protection_percentage, na.rm = TRUE), 1),
      median_protection_pct = round(median(protection_percentage, na.rm = TRUE), 1),
      min_protection_pct = round(min(protection_percentage, na.rm = TRUE), 1),
      max_protection_pct = round(max(protection_percentage, na.rm = TRUE), 1),
      sd_protection_pct = round(sd(protection_percentage, na.rm = TRUE), 1),
      mean_total_contribution = round(mean(total_contribution_percent, na.rm = TRUE), 1),
      n_years = n(),
      .groups = "drop"
    ) %>%
    arrange(mgmt_river)  # Maintains watershed ordering
  
  #--------------------------------------------------------------------------
  # EXPORT DATA
  #--------------------------------------------------------------------------
  
  if (export_csv) {
    export_closure_data(combined_results, summary_stats, output_dir, closure_start, closure_end)
  }
  
  cat("✓ Front-end closure protection analysis complete\n")
  cat(glue("  Output saved to: {output_dir}\n"))
  cat("  Files created:\n")
  cat("    - front_end_closure_protection_boxplot.png\n")
  if (export_csv) {
    cat(glue("    - closure_protection_{closure_start}_to_{closure_end}.csv\n"))
    cat(glue("    - closure_protection_summary_{closure_start}_to_{closure_end}.csv\n"))
    cat(glue("    - closure_protection_timeseries_{closure_start}_to_{closure_end}.csv\n"))
  }
  
  return(list(
    protection_data = combined_results,
    summary_stats = summary_stats,
    overall_protection = overall_protection
  ))
}

################################################################################
# DATA EXPORT FUNCTIONS
################################################################################

#' Export closure protection analysis data
export_closure_data <- function(closure_data, summary_stats, output_dir, closure_start, closure_end) {
  
  cat("Exporting closure protection data...\n")
  
  #--------------------------------------------------------------------------
  # EXPORT 1: DETAILED RESULTS
  #--------------------------------------------------------------------------
  
  detailed_results <- closure_data %>%
    select(
      year, mgmt_river, total_production, closure_production, 
      protection_percentage, total_contribution_percent, closure_window
    ) %>%
    arrange(year, desc(total_contribution_percent))
  
  detailed_filepath <- file.path(output_dir, glue("closure_protection_{closure_start}_to_{closure_end}.csv"))
  write_csv(detailed_results, detailed_filepath)
  cat(glue("Exported detailed results to: {basename(detailed_filepath)}\n"))
  
  #--------------------------------------------------------------------------
  # EXPORT 2: SUMMARY STATISTICS
  #--------------------------------------------------------------------------
  
  summary_filepath <- file.path(output_dir, glue("closure_protection_summary_{closure_start}_to_{closure_end}.csv"))
  write_csv(summary_stats, summary_filepath)
  cat(glue("Exported summary statistics to: {basename(summary_filepath)}\n"))
  
  #--------------------------------------------------------------------------
  # EXPORT 3: TIMESERIES FORMAT
  #--------------------------------------------------------------------------
  
  # Create timeseries for each management unit
  mgmt_timeseries <- NULL
  
  for (mgmt_name in unique(closure_data$mgmt_river)) {
    mgmt_data <- closure_data %>% filter(mgmt_river == mgmt_name)
    
    # Total contribution row
    total_row <- mgmt_data %>%
      select(year, total_contribution_percent) %>%
      pivot_wider(names_from = year, values_from = total_contribution_percent, values_fill = 0) %>%
      mutate(
        Management_Unit = mgmt_name,
        Metric = "Total_Contribution_Percent"
      )
    
    # Protection percentage row
    protection_row <- mgmt_data %>%
      select(year, protection_percentage) %>%
      pivot_wider(names_from = year, values_from = protection_percentage, values_fill = 0) %>%
      mutate(
        Management_Unit = mgmt_name,
        Metric = "Protection_Percentage"
      )
    
    # Combine both rows
    mgmt_timeseries <- bind_rows(mgmt_timeseries, total_row, protection_row)
  }
  
  # Sort by management unit, then metric
  mgmt_timeseries <- mgmt_timeseries %>%
    arrange(Management_Unit, Metric)
  
  timeseries_filepath <- file.path(output_dir, glue("closure_protection_timeseries_{closure_start}_to_{closure_end}.csv"))
  write_csv(mgmt_timeseries, timeseries_filepath)
  cat(glue("Exported timeseries format to: {basename(timeseries_filepath)}\n"))
  
  #--------------------------------------------------------------------------
  # EXPORT 4: OVERALL PROTECTION BY YEAR
  #--------------------------------------------------------------------------
  
  overall_protection <- closure_data %>%
    group_by(year) %>%
    summarise(
      total_annual_production = sum(total_production),
      total_closure_production = sum(closure_production),
      overall_protection_pct = (total_closure_production / total_annual_production) * 100,
      n_management_units = n(),
      .groups = "drop"
    )
  
  overall_filepath <- file.path(output_dir, glue("overall_protection_by_year_{closure_start}_to_{closure_end}.csv"))
  write_csv(overall_protection, overall_filepath)
  cat(glue("Exported overall protection by year to: {basename(overall_filepath)}\n"))
  
  cat("All closure protection CSV exports completed!\n")
}

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Quick closure protection summary
get_closure_summary <- function(mgmt_unit, year, closure_start = CONFIG$closure_dates$start, closure_end = CONFIG$closure_dates$end) {
  
  results_file <- file.path(PATHS$output_dir, "Closure_Protection", 
                            glue("closure_protection_{closure_start}_to_{closure_end}.csv"))
  
  if (file.exists(results_file)) {
    results <- read_csv(results_file, show_col_types = FALSE)
    unit_data <- results %>% filter(mgmt_river == mgmt_unit, year == !!year)
    
    if (nrow(unit_data) > 0) {
      cat(glue("Protection for {mgmt_unit} in {year}:\n"))
      cat(glue("  Protection percentage: {round(unit_data$protection_percentage, 1)}%\n"))
      return(unit_data)
    } else {
      cat(glue("No data found for {mgmt_unit} in {year}\n"))
      return(NULL)
    }
  } else {
    cat("No closure protection results found. Run run_closure_analysis() first.\n")
    return(NULL)
  }
}

################################################################################
# EXECUTION SECTION
################################################################################

# Run the analysis if this script is executed directly
if (interactive() || !exists(".closure_script_executed")) {
  cat("Front-end closure protection analysis functions loaded.\n")
  cat("Make sure you have run:\n")
  cat("  source('00_setup.R')\n")
  cat("  source('05_visualization_functions.R')\n")
  cat("before running this script.\n\n")
  
  cat("Available functions:\n")
  cat("  run_closure_analysis()                    # Main closure protection analysis\n")
  cat("  get_closure_summary(mgmt_unit, year)      # Quick summary for specific unit/year\n\n")
  
  # Uncomment the line below to run the main analysis
  # results <- run_closure_analysis()
  
  .closure_script_executed <- TRUE
}

cat("✓ Front-end closure protection analysis functions loaded\n")