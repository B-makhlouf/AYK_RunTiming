################################################################################
# MAIN.R - COMPLETE SALMON RUN TIMING ANALYSIS SUITE
################################################################################
# PURPOSE: Creates management river maps and exports CSV data:
#   1. DOY_Quartile: Shows proportion within each timing quartile
#   2. DOY_Total: Shows proportion of total annual run 
#   3. Front-end Closure: Shows protection during closure window
#   4. Average Quartile Maps: Shows average production across all years
#   5. Cumulative Distribution: Shows timing progression by management unit
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
library(scales)
library(readr)
library(png)

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
        HUC = 8,
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
        HUC = 8,
        return_values = FALSE  # We don't need return values for this
      )
    }
  }
  
  message("DOY Total Analysis complete - tributary maps created")
}

################################################################################
# ANALYSIS 3: FRONT-END CLOSURE BOXPLOT ANALYSIS
# PURPOSE: Shows what % of each management unit's production falls within Q1 closure
################################################################################

#' Run Front-end Closure Boxplot Analysis Only
run_closure_boxplot_analysis <- function(years, watershed = "Kusko") {
  
  message("=== Starting Front-End Closure Boxplot Analysis ===")
  
  # Output directory
  OUTPUT_DIR <- here("Figures/Front_End_Closure_Boxplots")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load ALL quartile data (Q1, Q2, Q3, Q4) to calculate totals
  DATA_PATH <- here("Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")
  
  if (!file.exists(DATA_PATH)) {
    warning("Management river analysis data not found. Run DOY analysis first.")
    return(NULL)
  }
  
  mgmt_data <- read.csv(DATA_PATH)
  
  # Remove Johnson river and clean quartile names
  mgmt_data_clean <- mgmt_data %>%
    filter(mgmt_river != "Johnson") %>%
    filter(year %in% years) %>%  # Filter to specified years
    mutate(quartile_clean = case_when(
      quartile == "Q1" ~ "Q1",
      quartile == "Q2" ~ "Q2", 
      quartile == "Q3" ~ "Q3",
      quartile == "Q4" ~ "Q4",
      TRUE ~ quartile
    ))
  
  message(paste("Loaded data for", length(unique(mgmt_data_clean$mgmt_river)), "management units"))
  message(paste("Years:", paste(sort(unique(mgmt_data_clean$year)), collapse = ", ")))
  
  # Calculate total annual production for each management unit in each year
  # Sum the total_run_prop across all quartiles (Q1+Q2+Q3+Q4)
  annual_totals <- mgmt_data_clean %>%
    group_by(year, mgmt_river) %>%
    summarise(
      total_annual_production = sum(total_run_prop, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Get Q1 production values
  q1_production <- mgmt_data_clean %>%
    filter(quartile_clean == "Q1") %>%
    select(year, mgmt_river, q1_production = total_run_prop)
  
  # Join Q1 with annual totals and calculate closure percentage
  boxplot_data <- q1_production %>%
    left_join(annual_totals, by = c("year", "mgmt_river")) %>%
    mutate(
      closure_percentage = (q1_production / total_annual_production) * 100
    ) %>%
    select(year, mgmt_river, q1_production, total_annual_production, closure_percentage)
  
  # Set watershed order (upstream to downstream)
  watershed_order <- c(
    "N. Fork Kusko", "E. Fork Kuskokwim", "S. Fork Kusko", 
    "Takotna and Nixon Fork", "Big River", "Upper Kusko Main",
    "Tatlawiksuk", "Kwethluk", "Stony", "Swift",
    "Holitna and Hoholitna", "George", "Oskakawlik", 
    "Middle Kusko Main", "Holokuk", "Aniak", "Tuluksak",
    "Kisaralik", "Hoholitna", "Johnson", "Lower Kusko"
  )
  
  # Filter to units present in data, maintain watershed order
  units_in_data <- unique(boxplot_data$mgmt_river)
  final_order <- watershed_order[watershed_order %in% units_in_data]
  
  # Apply ordering (reverse for coord_flip so upstream appears at top)
  boxplot_data$mgmt_river <- factor(boxplot_data$mgmt_river, levels = rev(final_order))
  
  year_range <- paste0(min(boxplot_data$year), "-", max(boxplot_data$year))
  
  # Create boxplot
  boxplot <- ggplot(boxplot_data, aes(x = mgmt_river, y = closure_percentage)) +
    geom_boxplot(
      fill = "lightblue", 
      alpha = 0.7, 
      outlier.size = 2.5,
      outlier.shape = 16,
      linewidth = 0.6,
      color = "darkblue"
    ) +
    coord_flip() +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, max(boxplot_data$closure_percentage, na.rm = TRUE) * 1.05),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      title = "Front-End Closure Protection by Management Unit",
      subtitle = paste("% of EACH UNIT'S total annual production within Q1 closure window | Years:", year_range),
      x = "Management Unit (Watershed Position: upstream → downstream)",
      y = "% of Unit's Total Annual Production in Q1 Closure Window",
      caption = "Shows what % of each unit's total annual production (Q1+Q2+Q3+Q4) occurs within Q1 closure period"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey20"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
  
  # Save boxplot
  ggsave(file.path(OUTPUT_DIR, "front_end_closure_protection_boxplot.png"), 
         boxplot, width = 12, height = 10, dpi = 300, bg = "white")
  
  # Summary statistics
  summary_stats <- boxplot_data %>%
    group_by(mgmt_river) %>%
    summarise(
      mean_closure_pct = round(mean(closure_percentage, na.rm = TRUE), 1),
      median_closure_pct = round(median(closure_percentage, na.rm = TRUE), 1),
      min_closure_pct = round(min(closure_percentage, na.rm = TRUE), 1),
      max_closure_pct = round(max(closure_percentage, na.rm = TRUE), 1),
      n_years = n(),
      .groups = "drop"
    ) %>%
    mutate(mgmt_river = factor(mgmt_river, levels = final_order)) %>%
    arrange(mgmt_river)
  
  write.csv(summary_stats, file.path(OUTPUT_DIR, "closure_protection_summary.csv"), row.names = FALSE)
  
  message(paste("Front-end closure boxplot completed. Saved to:", OUTPUT_DIR))
  message("Created:")
  message("  - front_end_closure_protection_boxplot.png")
  message("  - closure_protection_summary.csv")
  
  return(boxplot_data)
}

################################################################################
# ANALYSIS 4: AVERAGE QUARTILE MAPS
# PURPOSE: Shows average production across all years by quartile
################################################################################

#' Run Average Quartile Maps Analysis
run_average_quartile_analysis <- function(years, watershed = "Kusko") {
  
  message("=== Starting Average Quartile Maps Analysis ===")
  
  # Set file paths
  DATA_PATH <- here("Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")
  SPATIAL_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
  BASIN_PATH <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
  OUTPUT_DIR <- here("Figures/Average_Management_Production", paste0(min(years), "_", max(years)))
  
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load management unit production data
  if (!file.exists(DATA_PATH)) {
    warning("Management river analysis data not found. Run DOY analysis first.")
    return(NULL)
  }
  
  mgmt_data <- read_csv(DATA_PATH) %>%
    filter(mgmt_river != "Johnson") %>%  # Remove Johnson as in other analyses
    filter(year %in% years) %>%  # Filter to specified years
    mutate(
      quartile_clean = case_when(
        quartile == "Q1" ~ "Q1",
        quartile == "Q2" ~ "Q2", 
        quartile == "Q3" ~ "Q3",
        quartile == "Q4" ~ "Q4",
        TRUE ~ quartile
      )
    )
  
  # Check what years are available
  available_years <- sort(unique(mgmt_data$year))
  year_range <- paste0(min(available_years), "-", max(available_years))
  
  message(paste("Creating average maps for years:", year_range))
  message(paste("Found", length(unique(mgmt_data$mgmt_river)), "management units"))
  
  # Load spatial data
  edges <- st_read(SPATIAL_PATH, quiet = TRUE)
  basin <- st_read(BASIN_PATH, quiet = TRUE)
  
  # Calculate average production by quartile
  avg_production_by_quartile <- mgmt_data %>%
    group_by(mgmt_river, quartile_clean) %>%
    summarise(
      avg_within_quartile_prop = mean(within_quartile_prop, na.rm = TRUE),
      n_years = n(),
      .groups = "drop"
    ) %>%
    # Convert to the same scale as original maps (0-1 proportion)
    mutate(production_proportion = avg_within_quartile_prop)
  
  # Create maps for each quartile
  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    message(paste("Creating average map for", q))
    
    # Get average data for this quartile
    quartile_data <- avg_production_by_quartile %>%
      filter(quartile_clean == q) %>%
      filter(!is.na(production_proportion))
    
    # Create map with exact original styling
    map_filename <- paste0("Average_", q, "_Management_", year_range, ".png")
    map_filepath <- file.path(OUTPUT_DIR, map_filename)
    
    create_average_mgmt_map(quartile_data, q, map_filepath, year_range, edges, basin)
  }
  
  # Save average production data
  write_csv(avg_production_by_quartile, 
            file.path(OUTPUT_DIR, paste0("average_production_by_quartile_", year_range, ".csv")))
  
  message(paste("Average quartile maps completed. Saved to:", OUTPUT_DIR))
}

#' Create management map identical to original DOY_Quartile/Management style
create_average_mgmt_map <- function(quartile_data, quartile_label, output_filepath, year_range, edges, basin) {
  
  # Create the PNG with dimensions for single map (no bar plot)
  png(file = output_filepath, width = 10, height = 8, units = "in", res = 300, bg = "white")
  
  # Join with spatial data and set up linewidths like DFA maps
  edges_with_data <- edges %>%
    left_join(quartile_data, by = "mgmt_river") %>%
    filter(!is.na(mgmt_river) & mgmt_river != "") %>%
    mutate(
      # Set production proportion (handle NA values)
      production_proportion = ifelse(is.na(production_proportion), 0, production_proportion),
      # Set stream order and linewidth exactly like DFA maps
      stream_order = ifelse(is.na(Str_Order), 3, Str_Order),
      line_width = pmax(0.3, pmin(3.0, 0.3 + (stream_order - min(stream_order, na.rm = TRUE)) * 
                                    (3.0 - 0.3) / (max(stream_order, na.rm = TRUE) - min(stream_order, na.rm = TRUE))))
    )
  
  # Ensure consistent CRS
  if (st_crs(basin) != st_crs(edges_with_data)) {
    basin <- st_transform(basin, st_crs(edges_with_data))
  }
  
  # Main map plot - single map without bar chart
  main_plot <- ggplot() +
    geom_sf(data = basin, fill = "gray95", color = "gray70", 
            linewidth = 0.5, alpha = 0.3) +
    geom_sf(data = edges_with_data, 
            aes(color = production_proportion, linewidth = stream_order), 
            alpha = 0.8) +
    scale_color_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Production\nProportion",
      na.value = "grey60",
      labels = scales::percent_format(accuracy = 1),
      guide = guide_colorbar(
        barwidth = 1, barheight = 15,
        frame.colour = "grey40", ticks.colour = "grey40",
        show.limits = TRUE
      )
    ) +
    scale_linewidth_continuous(
      range = c(0.3, 3.0), 
      name = "Stream\nOrder"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Average ", quartile_label, ": Management Rivers - Kusko Watershed"),
      subtitle = paste("Average across all years (", year_range, ")", sep = "")
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "grey30"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey50"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold", color = "grey30"),
      legend.text = element_text(color = "grey30"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10, "mm")
    )
  
  print(main_plot)
  
  dev.off()
  
  message(paste("Created average map:", basename(output_filepath)))
}

################################################################################
# ANALYSIS 5: CUMULATIVE DISTRIBUTION ANALYSIS
# PURPOSE: Shows timing progression by management unit with 3-day intervals
################################################################################

#' Run Cumulative Distribution Analysis
run_cumulative_distribution_analysis <- function(years, watershed = "Kusko", interval_days = 3) {
  
  message("=== Starting Cumulative Distribution Analysis ===")
  
  # Source the cumulative distribution functions
  if (file.exists(here("Code/analysis/Cumulative_Distribution.R"))) {
    # Extract just the function definitions from the cumulative distribution script
    source(here("Code/analysis/Cumulative_Distribution.R"))
    
    # Run the analysis
    results <- analyze_cumulative_distributions(
      years = years,
      watershed = watershed,
      interval_days = interval_days,
      export_csv = TRUE
    )
    
    if (!is.null(results)) {
      message("Cumulative distribution analysis completed successfully")
      return(results)
    } else {
      warning("Cumulative distribution analysis failed")
      return(NULL)
    }
  } else {
    warning("Cumulative_Distribution.R not found. Skipping cumulative analysis.")
    return(NULL)
  }
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
years <- c(2017, 2018, 2019, 2020, 2021)
watersheds <- c("Kusko")

message("=== Running Complete Analysis Suite ===")

# 1. DOY Quartiles (Management Rivers only) + CSV export
message("1. Starting DOY Quartile Analysis (Management Rivers only) with CSV export...")
run_doy_analysis(years, watersheds, export_csv = TRUE)

# # 2. DOY Total (tributary maps only - no HUC maps)
# message("2. Starting DOY Total Analysis (Tributary maps only)...")
# run_doy_total_analysis(years, watersheds, export_csv = FALSE)

# 3. Front-end Closure Boxplot Analysis
message("3. Starting Front-end Closure Boxplot Analysis...")
run_closure_boxplot_analysis(years, "Kusko")

# 4. Average Quartile Maps
message("4. Starting Average Quartile Maps Analysis...")
run_average_quartile_analysis(years, "Kusko")

# 5. Cumulative Distribution Analysis
message("5. Starting Cumulative Distribution Analysis...")
run_cumulative_distribution_analysis(years, "Kusko", interval_days = 3)

message("=== Complete Analysis Suite Finished ===")
message("All analyses have been completed!")
message("Check the following directories for outputs:")
message("  - Basin Maps/DOY_Quartile/")
message("  - Basin Maps/DOY_Total/")
message("  - Figures/Front_End_Closure_Boxplots/")
message("  - Figures/Average_Management_Production/")
message("  - Analysis_Results/Cumulative_Distribution/")
message("  - Analysis_Results/Management_River_Analysis/")