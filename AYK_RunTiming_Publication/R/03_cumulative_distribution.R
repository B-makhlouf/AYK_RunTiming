################################################################################
# 03_CUMULATIVE_DISTRIBUTION.R - CUMULATIVE TIMING DISTRIBUTION ANALYSIS
################################################################################
# Shows how production builds up over time for each management unit
# Creates timing progression plots and run duration analysis
# Run 00_setup.R and 05_visualization_functions.R first
################################################################################

cat("=== CUMULATIVE DISTRIBUTION ANALYSIS ===\n")

# Check if setup is loaded
if (!exists("CONFIG")) {
  stop("Please run 00_setup.R first to load parameters and functions")
}

if (!exists("create_cumulative_plots")) {
  stop("Please run 05_visualization_functions.R first to load plotting functions")
}

################################################################################
# MAIN CUMULATIVE DISTRIBUTION ANALYSIS FUNCTION
################################################################################

#' Run cumulative distribution analysis
#' Shows how production builds up over time for each management unit
#' Produces identical outputs to original cumulative distribution analysis
run_cumulative_analysis <- function(years = CONFIG$years,
                                    watersheds = CONFIG$watersheds,
                                    interval_days = CONFIG$cumulative_interval_days,
                                    export_csv = TRUE) {
  
  cat("Starting cumulative distribution analysis...\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n")
  cat("Data points every", interval_days, "days\n")
  
  # Create output directory
  output_dir <- file.path(PATHS$output_dir, "Cumulative_Distribution")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_cumulative_data <- list()
  all_run_duration_data <- list()
  
  for (watershed in watersheds) {
    params <- WATERSHED_PARAMS[[watershed]]
    
    for (year in years) {
      cat(glue("Processing cumulative distribution for {watershed} {year}...\n"))
      
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
      
      # Create DOY intervals every N days
      doy_range <- range(natal_data$DOY, na.rm = TRUE)
      doy_sequence <- seq(from = floor(doy_range[1]), 
                          to = ceiling(doy_range[2]), 
                          by = interval_days)
      
      # Ensure we include the final DOY
      if (max(doy_sequence) < ceiling(doy_range[2])) {
        doy_sequence <- c(doy_sequence, ceiling(doy_range[2]))
      }
      
      cat(glue("  DOY range: {floor(doy_range[1])} to {ceiling(doy_range[2])}\n"))
      cat(glue("  Created {length(doy_sequence)} intervals every {interval_days} days\n"))
      
      year_cumulative_data <- data.frame()
      
      # Process each time interval
      for (i in 1:length(doy_sequence)) {
        current_doy <- doy_sequence[i]
        
        # Get cumulative data up to this DOY
        cumulative_subset <- natal_data %>% filter(DOY <= current_doy)
        
        if (nrow(cumulative_subset) == 0) next
        
        cat(glue("    Processing DOY {current_doy} with {nrow(cumulative_subset)} fish\n"))
        
        # Perform assignment
        assignment_matrix <- perform_assignment(
          cumulative_subset, spatial_data$edges, watershed, priors,
          pid_iso, error, params$sensitivity_threshold
        )
        
        # Process by management unit
        basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
        mgmt_result <- process_mgmt_assignments(spatial_data$edges, basin_assign_sum)
        
        if (!is.null(mgmt_result)) {
          # Add temporal information
          interval_data <- mgmt_result %>%
            mutate(
              year = year,
              doy = current_doy,
              interval_number = i,
              cumulative_production = total_production,
              fish_count = nrow(cumulative_subset)
            )
          
          year_cumulative_data <- bind_rows(year_cumulative_data, interval_data)
        }
      }
      
      # Calculate percentages within each management unit
      if (nrow(year_cumulative_data) > 0) {
        year_cumulative_data <- year_cumulative_data %>%
          group_by(mgmt_river) %>%
          arrange(doy) %>%
          mutate(
            final_production = max(cumulative_production),
            cumulative_percent = (cumulative_production / final_production) * 100,
            total_production_all_units = sum(final_production)
          ) %>%
          ungroup() %>%
          mutate(
            overall_percent = (cumulative_production / total_production_all_units) * 100
          )
        
        # Calculate run duration metrics
        year_run_duration <- year_cumulative_data %>%
          group_by(mgmt_river) %>%
          arrange(doy) %>%
          summarise(
            year = first(year),
            run_start_doy = min(doy, na.rm = TRUE),
            
            # Find first DOY where cumulative percent goes above 50%
            doy_above_50_percent = {
              above_50_idx <- which(cumulative_percent > 50)
              if (length(above_50_idx) > 0) {
                doy[min(above_50_idx)]
              } else {
                max(doy, na.rm = TRUE)
              }
            },
            
            # Actual percentage at that DOY
            actual_percent_at_doy = {
              above_50_idx <- which(cumulative_percent > 50)
              if (length(above_50_idx) > 0) {
                cumulative_percent[min(above_50_idx)]
              } else {
                max(cumulative_percent, na.rm = TRUE)
              }
            },
            
            # Days from run start to above 50%
            days_to_50_percent = doy_above_50_percent - run_start_doy,
            
            # Additional metrics
            max_percent_reached = max(cumulative_percent, na.rm = TRUE),
            total_production = max(cumulative_production, na.rm = TRUE),
            final_overall_percent = max(overall_percent, na.rm = TRUE),
            .groups = "drop"
          )
        
        all_cumulative_data[[glue("{watershed}_{year}")]] <- year_cumulative_data
        all_run_duration_data[[glue("{watershed}_{year}")]] <- year_run_duration
      }
      
      cat(glue("  Completed year {year} with {nrow(year_cumulative_data)} data points\n"))
    }
  }
  
  #--------------------------------------------------------------------------
  # COMBINE ALL YEARS AND CREATE VISUALIZATIONS
  #--------------------------------------------------------------------------
  
  if (length(all_cumulative_data) == 0) {
    stop("No valid data processed for any year")
  }
  
  # Combine all data
  combined_cumulative_data <- bind_rows(all_cumulative_data)
  combined_run_duration <- bind_rows(all_run_duration_data)
  
  cat(glue("Combined data: {nrow(combined_cumulative_data)} total data points\n"))
  cat(glue("Management units analyzed: {length(unique(combined_cumulative_data$mgmt_river))}\n"))
  
  # Print management units found
  final_mgmt_units <- sort(unique(combined_cumulative_data$mgmt_river))
  cat("\nManagement units in final data:\n")
  for (i in 1:length(final_mgmt_units)) {
    cat(glue("  {i}. '{final_mgmt_units[i]}'\n"))
  }
  
  #--------------------------------------------------------------------------
  # CREATE VISUALIZATION PLOTS
  #--------------------------------------------------------------------------
  
  create_cumulative_plots(combined_cumulative_data, output_dir, watersheds[1], interval_days)
  
  #--------------------------------------------------------------------------
  # CREATE RUN DURATION ANALYSIS
  #--------------------------------------------------------------------------
  
  cat("Creating run duration analysis...\n")
  
  if (nrow(combined_run_duration) > 0) {
    # Apply watershed ordering
    combined_run_duration <- apply_watershed_order(combined_run_duration, "mgmt_river", reverse_for_plots = TRUE)
    
    # [trimmed: not a paper figure] run-duration range plot not exported.
    # create_run_duration_range_plot(combined_run_duration, output_dir, watersheds[1], interval_days)
    
    # Calculate and display summary statistics
    duration_summary <- combined_run_duration %>%
      group_by(mgmt_river) %>%
      summarise(
        avg_days_to_50 = round(mean(days_to_50_percent, na.rm = TRUE), 1),
        min_days = min(days_to_50_percent, na.rm = TRUE),
        max_days = max(days_to_50_percent, na.rm = TRUE),
        avg_actual_percent = round(mean(actual_percent_at_doy, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>%
      arrange(avg_days_to_50)
    
    cat("\nRun Duration Summary (days from start to 50% completion):\n")
    print(duration_summary)
    
    if (nrow(duration_summary) > 0) {
      fastest_unit <- duration_summary$mgmt_river[1]
      slowest_unit <- duration_summary$mgmt_river[nrow(duration_summary)]
      fastest_days <- duration_summary$avg_days_to_50[1]
      slowest_days <- duration_summary$avg_days_to_50[nrow(duration_summary)]
      
      cat(glue("\n🏃 FASTEST BUILD-UP: {fastest_unit} - avg {fastest_days} days to 50%\n"))
      cat(glue("🐌 SLOWEST BUILD-UP: {slowest_unit} - avg {slowest_days} days to 50%\n"))
    }
  }
  
  #--------------------------------------------------------------------------
  # EXPORT DATA
  #--------------------------------------------------------------------------
  
  if (export_csv) {
    export_cumulative_data(combined_cumulative_data, combined_run_duration, 
                           output_dir, watersheds[1], interval_days)
  }
  
  cat("✓ Cumulative distribution analysis complete\n")
  cat(glue("  Output saved to: {output_dir}\n"))
  
  return(list(
    cumulative_data = combined_cumulative_data,
    run_duration_data = combined_run_duration,
    parameters = list(
      watershed = watersheds[1],
      years = years,
      interval_days = interval_days
    )
  ))
}

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Create run duration range plot (horizontal orientation)
create_run_duration_range_plot <- function(run_duration_data, output_dir, watershed, interval_days) {
  
  plot_dir <- file.path(output_dir, "Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Filter clean data
  clean_data <- run_duration_data %>%
    filter(!is.na(mgmt_river), !is.na(days_to_50_percent))
  
  if (nrow(clean_data) == 0) {
    cat("❌ No clean data available for range plot\n")
    return(NULL)
  }
  
  # Get management units in watershed order
  mgmt_rivers <- sort(unique(clean_data$mgmt_river))
  final_order <- WATERSHED_ORDER[WATERSHED_ORDER %in% mgmt_rivers]
  missing_units <- setdiff(mgmt_rivers, WATERSHED_ORDER)
  if (length(missing_units) > 0) {
    final_order <- c(final_order, missing_units)
  }
  
  clean_data$mgmt_river <- factor(clean_data$mgmt_river, levels = final_order)
  
  # Create colors matching cumulative plots
  n_units <- length(final_order)
  watershed_colors <- colorRampPalette(c(
    "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
    "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
  ))(n_units)
  names(watershed_colors) <- final_order
  
  # Calculate summary statistics
  duration_stats <- clean_data %>%
    group_by(mgmt_river) %>%
    summarise(
      n_years = n(),
      mean_days = mean(days_to_50_percent, na.rm = TRUE),
      min_days = min(days_to_50_percent, na.rm = TRUE),
      max_days = max(days_to_50_percent, na.rm = TRUE),
      range_days = max_days - min_days,
      .groups = 'drop'
    ) %>%
    arrange(mean_days)
  
  duration_stats$mgmt_river <- factor(duration_stats$mgmt_river, levels = final_order)
  
  # Create modern horizontal range plot
  p_range <- ggplot(duration_stats, aes(y = mgmt_river)) +
    # Range lines (min to max)
    geom_segment(aes(x = min_days, xend = max_days, color = mgmt_river), 
                 linewidth = 6, alpha = 0.7) +
    
    # Mean markers (black vertical lines)
    geom_point(aes(x = mean_days), color = "black", size = 4, shape = "|", stroke = 2) +
    
    # Colors matching cumulative plots
    scale_color_manual(values = watershed_colors, breaks = final_order, guide = "none") +
    
    scale_x_continuous(
      name = "Days to 50% return",
      breaks = function(x) pretty(x, n = 6),
      expand = expansion(mult = c(0.02, 0.02)),
      labels = function(x) paste0(x, " days")
    ) +
    
    scale_y_discrete(name = NULL, limits = rev(final_order)) +
    
    labs(title = "Days to 50% of return by stock") +
    
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0, color = "gray15", margin = margin(b = 5)),
      axis.text.y = element_text(size = 11, color = "gray30", hjust = 1, margin = margin(r = 10)),
      axis.text.x = element_text(size = 10, color = "gray30", margin = margin(t = 8)),
      axis.title.x = element_text(size = 12, color = "gray20", face = "bold", margin = margin(t = 15)),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.4, linetype = "solid"),
      panel.grid.minor.x = element_line(color = "gray96", linewidth = 0.2),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(20, 25, 20, 20, "mm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Save plot
  range_filename <- glue("run_duration_range_{watershed}_{interval_days}day_intervals.png")
  ggsave(file.path(plot_dir, range_filename), p_range, 
         width = 12, height = 8, dpi = 300, bg = "white")
  
  cat(glue("✅ Created run duration range plot: {range_filename}\n"))
  
  return(duration_stats)
}

#' Export cumulative distribution data
export_cumulative_data <- function(cumulative_data, run_duration_data, output_dir, watershed, interval_days) {
  
  csv_dir <- file.path(output_dir, "CSV")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Exporting cumulative distribution data...\n")
  
  # Export full cumulative data
  cumulative_clean <- cumulative_data %>%
    select(year, doy, interval_number, mgmt_river, cumulative_production, 
           cumulative_percent, overall_percent, fish_count) %>%
    arrange(year, mgmt_river, doy)
  
  cumulative_filepath <- file.path(csv_dir, glue("cumulative_distribution_data_{watershed}_{interval_days}day_intervals.csv"))
  write_csv(cumulative_clean, cumulative_filepath)
  cat(glue("Exported cumulative data to: {basename(cumulative_filepath)}\n"))
  
  # Export run duration metrics
  run_duration_clean <- run_duration_data %>%
    select(year, mgmt_river, run_start_doy, doy_above_50_percent, actual_percent_at_doy, 
           days_to_50_percent, max_percent_reached, total_production, final_overall_percent) %>%
    arrange(year, mgmt_river)
  
  run_duration_filepath <- file.path(csv_dir, glue("run_duration_metrics_{watershed}_{interval_days}day_intervals.csv"))
  write_csv(run_duration_clean, run_duration_filepath)
  cat(glue("Exported run duration metrics to: {basename(run_duration_filepath)}\n"))
  
  # Export wide format for easy plotting
  cumulative_wide <- cumulative_data %>%
    select(year, doy, mgmt_river, cumulative_percent) %>%
    pivot_wider(names_from = mgmt_river, values_from = cumulative_percent, 
                values_fill = 0, names_prefix = "mgmt_") %>%
    arrange(year, doy)
  
  wide_filepath <- file.path(csv_dir, glue("cumulative_distribution_wide_{watershed}_{interval_days}day_intervals.csv"))
  write_csv(cumulative_wide, wide_filepath)
  cat(glue("Exported wide format data to: {basename(wide_filepath)}\n"))
  
  # Export summary statistics
  summary_stats <- cumulative_data %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = max(cumulative_production, na.rm = TRUE),
      final_overall_percent = max(overall_percent, na.rm = TRUE),
      n_data_points = n(),
      first_doy = min(doy, na.rm = TRUE),
      last_doy = max(doy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_production))
  
  summary_filepath <- file.path(csv_dir, glue("timing_summary_stats_{watershed}_{interval_days}day_intervals.csv"))
  write_csv(summary_stats, summary_filepath)
  cat(glue("Exported summary statistics to: {basename(summary_filepath)}\n"))
  
  cat("All CSV exports completed!\n")
}

################################################################################
# SPECIALIZED ANALYSIS FUNCTIONS
################################################################################

#' Run just the timing progression analysis (no run duration)
run_timing_progression_only <- function(years = CONFIG$years, watersheds = CONFIG$watersheds, interval_days = CONFIG$cumulative_interval_days) {
  
  cat("Running timing progression analysis only...\n")
  
  # Run full analysis but only return cumulative data
  results <- run_cumulative_analysis(years, watersheds, interval_days, export_csv = FALSE)
  
  # Create only the timing plots
  output_dir <- file.path(PATHS$output_dir, "Timing_Progression_Only")
  create_cumulative_plots(results$cumulative_data, output_dir, watersheds[1], interval_days)
  
  cat(glue("✓ Timing progression analysis complete. Output in: {output_dir}\n"))
  
  return(results$cumulative_data)
}

#' Run just the run duration analysis (requires existing cumulative data)
run_duration_analysis_only <- function(years = CONFIG$years, watersheds = CONFIG$watersheds, interval_days = CONFIG$cumulative_interval_days) {
  
  cat("Running run duration analysis only...\n")
  
  # Check if cumulative data exists
  cumulative_file <- file.path(PATHS$output_dir, "Cumulative_Distribution/CSV", 
                               glue("run_duration_metrics_{watersheds[1]}_{interval_days}day_intervals.csv"))
  
  if (file.exists(cumulative_file)) {
    # Load existing data
    run_duration_data <- read_csv(cumulative_file, show_col_types = FALSE)
    run_duration_data <- apply_watershed_order(run_duration_data, "mgmt_river", reverse_for_plots = TRUE)
    
    output_dir <- file.path(PATHS$output_dir, "Run_Duration_Only")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create duration plots
    duration_stats <- create_run_duration_range_plot(run_duration_data, output_dir, watersheds[1], interval_days)
    
    cat(glue("✓ Run duration analysis complete. Output in: {output_dir}\n"))
    return(duration_stats)
    
  } else {
    cat("No existing run duration data found. Running full cumulative analysis...\n")
    results <- run_cumulative_analysis(years, watersheds, interval_days, export_csv = TRUE)
    return(results$run_duration_data)
  }
}

################################################################################
# EXECUTION SECTION
################################################################################

# Run the analysis if this script is executed directly
if (interactive() || !exists(".cumulative_script_executed")) {
  cat("Cumulative distribution analysis functions loaded.\n")
  cat("Make sure you have run:\n")
  cat("  source('00_setup.R')\n")
  cat("  source('05_visualization_functions.R')\n")
  cat("before running this script.\n\n")
  
  cat("Available functions:\n")
  cat("  run_cumulative_analysis()        # Full cumulative distribution analysis\n")
  cat("  run_timing_progression_only()    # Just timing progression plots\n")
  cat("  run_duration_analysis_only()     # Just run duration analysis\n\n")
  
  # Uncomment the line below to run the full analysis
  # results <- run_cumulative_analysis()
  
  .cumulative_script_executed <- TRUE
}

cat("✓ Cumulative distribution analysis functions loaded\n")