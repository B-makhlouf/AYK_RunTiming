################################################################################
# 02_AVERAGE_MAPS.R - AVERAGE PRODUCTION MAPS ANALYSIS
################################################################################
# Creates average production maps by quartile & management unit with variability boxplots
# Requires 01_tributary_maps.R to be run first to generate underlying data
# Run 00_setup.R and 05_visualization_functions.R first
################################################################################

cat("=== AVERAGE PRODUCTION MAPPING ANALYSIS ===\n")

# Check if setup is loaded
if (!exists("CONFIG")) {
  stop("Please run 00_setup.R first to load parameters and functions")
}

if (!exists("create_average_production_map")) {
  stop("Please run 05_visualization_functions.R first to load plotting functions")
}

################################################################################
# MAIN AVERAGE PRODUCTION ANALYSIS FUNCTION
################################################################################

#' Run average production mapping analysis
#' Creates maps showing average production patterns across years with variability boxplots
#' Produces identical outputs to original average production analysis
run_average_maps_analysis <- function(years = CONFIG$years,
                                      watersheds = CONFIG$watersheds) {
  
  cat("Starting average production mapping analysis...\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n")
  
  for (watershed in watersheds) {
    year_range <- paste0(min(years), "-", max(years))
    
    cat(glue("Processing {watershed} for years {year_range}...\n"))
    
    # Check for management river data
    mgmt_data_path <- file.path(PATHS$output_dir, "Management_River_Analysis/management_river_analysis_tidy.csv")
    
    if (!file.exists(mgmt_data_path)) {
      cat("Management river analysis data not found. Running tributary analysis first...\n")
      
      # Check if tributary analysis function is available
      if (!exists("run_tributary_analysis")) {
        stop("Please run 01_tributary_maps.R first to generate the required data")
      }
      
      # Run tributary analysis to generate data
      run_tributary_analysis(years, watersheds, export_csv = TRUE)
    }
    
    # Load management unit production data
    mgmt_data <- read_csv(mgmt_data_path, show_col_types = FALSE) %>%
      filter(mgmt_river != "Johnson") %>%  # Remove Johnson as in original
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
    
    cat(glue("Loaded data for {length(unique(mgmt_data$mgmt_river))} management units\n"))
    
    # Load spatial data
    spatial_data <- load_spatial_data(watershed)
    
    # Create output directory
    output_dir <- file.path(PATHS$figures_dir, "Average_Management_Production", year_range)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    #--------------------------------------------------------------------------
    # PART 1: CREATE AVERAGE MAPS
    #--------------------------------------------------------------------------
    
    cat("Creating average production maps...\n")
    
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
    
    # Apply watershed ordering
    avg_production_by_quartile <- apply_watershed_order(avg_production_by_quartile, "mgmt_river", reverse_for_plots = FALSE)
    
    # Create maps for each quartile
    for (q in CONFIG$quartiles) {
      cat(glue("  Creating average map for {q}\n"))
      
      # Get average data for this quartile
      quartile_data <- avg_production_by_quartile %>%
        filter(quartile_clean == q) %>%
        filter(!is.na(production_proportion))
      
      # Create map with exact original styling
      map_filename <- glue("Average_{q}_Management_{year_range}.png")
      map_path <- file.path(output_dir, map_filename)
      
      create_average_production_map(
        quartile_data, q, spatial_data$edges, spatial_data$basin, 
        year_range, map_path
      )
    }
    
    # Save average production data
    write_csv(avg_production_by_quartile, 
              file.path(output_dir, glue("average_production_by_quartile_{year_range}.csv")))
    
    #--------------------------------------------------------------------------
    # PART 2: CREATE VARIABILITY BOXPLOTS
    #--------------------------------------------------------------------------
    
    cat("Creating variability boxplots...\n")
    
    create_production_boxplots(mgmt_data, watershed, years, output_dir)
    
    #--------------------------------------------------------------------------
    # PART 3: SAVE SUMMARY STATISTICS
    #--------------------------------------------------------------------------
    
    cat("Calculating summary statistics...\n")
    
    # Prepare data for statistics
    quartile_boxplot_data <- mgmt_data %>%
      select(year, mgmt_river, quartile_clean, within_quartile_prop) %>%
      mutate(within_quartile_pct = within_quartile_prop * 100)
    
    # Apply watershed ordering
    quartile_boxplot_data <- apply_watershed_order(quartile_boxplot_data, "mgmt_river", reverse_for_plots = FALSE)
    
    # Calculate summary statistics
    quartile_summary_stats <- quartile_boxplot_data %>%
      group_by(mgmt_river, quartile_clean) %>%
      summarise(
        mean_quartile_production = mean(within_quartile_pct, na.rm = TRUE),
        median_quartile_production = median(within_quartile_pct, na.rm = TRUE),
        sd_quartile_production = sd(within_quartile_pct, na.rm = TRUE),
        min_quartile_production = min(within_quartile_pct, na.rm = TRUE),
        max_quartile_production = max(within_quartile_pct, na.rm = TRUE),
        cv_quartile_production = sd_quartile_production / mean_quartile_production,
        n_years = n(),
        .groups = "drop"
      ) %>%
      arrange(mgmt_river, quartile_clean)
    
    write_csv(quartile_summary_stats, 
              file.path(output_dir, "management_unit_quartile_summary_statistics.csv"))
    
    # Calculate detailed statistics
    detailed_stats <- quartile_summary_stats %>%
      rename(
        mean_pct = mean_quartile_production,
        median_pct = median_quartile_production,
        sd_pct = sd_quartile_production,
        min_pct = min_quartile_production,
        max_pct = max_quartile_production,
        cv_pct = cv_quartile_production
      ) %>%
      arrange(quartile_clean, mgmt_river)
    
    write_csv(detailed_stats, 
              file.path(output_dir, "management_unit_detailed_statistics.csv"))
    
    cat("✓ Average production mapping analysis complete\n")
    cat("  Output saved to:\n")
    cat(glue("    {output_dir}\n"))
    cat("  Files created:\n")
    cat("    MAPS:\n")
    for (q in CONFIG$quartiles) {
      cat(glue("      - Average_{q}_Management_{year_range}.png\n"))
    }
    cat("    BOXPLOTS:\n")
    for (q in CONFIG$quartiles) {
      cat(glue("      - management_unit_{q}_boxplot.png\n"))
    }
    cat("      - management_unit_all_quartiles_boxplot.png\n")
    cat("    DATA:\n")
    cat(glue("      - average_production_by_quartile_{year_range}.csv\n"))
    cat("      - management_unit_quartile_summary_statistics.csv\n")
    cat("      - management_unit_detailed_statistics.csv\n")
  }
}

################################################################################
# HELPER FUNCTIONS FOR SPECIFIC OUTPUTS
################################################################################

#' Create just the average maps (no boxplots)
create_average_maps_only <- function(years = CONFIG$years, watersheds = CONFIG$watersheds) {
  
  cat("Creating average maps only (no boxplots)...\n")
  
  for (watershed in watersheds) {
    year_range <- paste0(min(years), "-", max(years))
    
    # Load data
    mgmt_data_path <- file.path(PATHS$output_dir, "Management_River_Analysis/management_river_analysis_tidy.csv")
    
    if (!file.exists(mgmt_data_path)) {
      stop("Management river analysis data not found. Run 01_tributary_maps.R first.")
    }
    
    mgmt_data <- read_csv(mgmt_data_path, show_col_types = FALSE) %>%
      filter(mgmt_river != "Johnson", year %in% years) %>%
      mutate(quartile_clean = case_when(
        quartile == "Q1" ~ "Q1", quartile == "Q2" ~ "Q2", 
        quartile == "Q3" ~ "Q3", quartile == "Q4" ~ "Q4", TRUE ~ quartile
      ))
    
    spatial_data <- load_spatial_data(watershed)
    
    # Calculate averages
    avg_production_by_quartile <- mgmt_data %>%
      group_by(mgmt_river, quartile_clean) %>%
      summarise(avg_within_quartile_prop = mean(within_quartile_prop, na.rm = TRUE), 
                .groups = "drop") %>%
      mutate(production_proportion = avg_within_quartile_prop)
    
    avg_production_by_quartile <- apply_watershed_order(avg_production_by_quartile, "mgmt_river")
    
    # Create output directory
    output_dir <- file.path(PATHS$figures_dir, "Average_Maps_Only", year_range)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create maps
    for (q in CONFIG$quartiles) {
      quartile_data <- avg_production_by_quartile %>%
        filter(quartile_clean == q, !is.na(production_proportion))
      
      map_filename <- glue("Average_{q}_Management_{year_range}.png")
      map_path <- file.path(output_dir, map_filename)
      
      create_average_production_map(
        quartile_data, q, spatial_data$edges, spatial_data$basin, 
        year_range, map_path
      )
    }
    
    cat(glue("✓ Average maps created for {watershed} in {output_dir}\n"))
  }
}

#' Create just the boxplots (no maps)
create_boxplots_only <- function(years = CONFIG$years, watersheds = CONFIG$watersheds) {
  
  cat("Creating boxplots only (no maps)...\n")
  
  for (watershed in watersheds) {
    year_range <- paste0(min(years), "-", max(years))
    
    # Load data
    mgmt_data_path <- file.path(PATHS$output_dir, "Management_River_Analysis/management_river_analysis_tidy.csv")
    
    if (!file.exists(mgmt_data_path)) {
      stop("Management river analysis data not found. Run 01_tributary_maps.R first.")
    }
    
    mgmt_data <- read_csv(mgmt_data_path, show_col_types = FALSE) %>%
      filter(mgmt_river != "Johnson", year %in% years) %>%
      mutate(quartile_clean = case_when(
        quartile == "Q1" ~ "Q1", quartile == "Q2" ~ "Q2", 
        quartile == "Q3" ~ "Q3", quartile == "Q4" ~ "Q4", TRUE ~ quartile
      ))
    
    # Create output directory
    output_dir <- file.path(PATHS$figures_dir, "Boxplots_Only", year_range)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create boxplots
    create_production_boxplots(mgmt_data, watershed, years, output_dir)
    
    cat(glue("✓ Boxplots created for {watershed} in {output_dir}\n"))
  }
}

################################################################################
# EXECUTION SECTION
################################################################################

# Run the analysis if this script is executed directly
if (interactive() || !exists(".average_script_executed")) {
  cat("Average production mapping functions loaded.\n")
  cat("Make sure you have run:\n")
  cat("  source('00_setup.R')\n")
  cat("  source('05_visualization_functions.R')\n")
  cat("  source('01_tributary_maps.R')  # Or run run_tributary_analysis() first\n")
  cat("before running this script.\n\n")
  
  cat("Available functions:\n")
  cat("  run_average_maps_analysis()  # Full analysis with maps and boxplots\n")
  cat("  create_average_maps_only()   # Just the average maps\n")
  cat("  create_boxplots_only()       # Just the variability boxplots\n\n")
  
  # Uncomment the line below to run the full analysis
  # run_average_maps_analysis()
  
  .average_script_executed <- TRUE
}

cat("✓ Average production mapping functions loaded\n")