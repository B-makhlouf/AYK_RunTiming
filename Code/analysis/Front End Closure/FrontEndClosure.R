################################################################################
#                     ORIGINAL MAP CREATION FUNCTIONS                         #
################################################################################

create_mgmt_closure_map <- function(window_result, total_result, edges, basin, year, watershed, 
                                    start_date, end_date, overall_pct, map_output_dir, 
                                    natal_data, window_data, total_production_all) {
  
  start_doy <- as.numeric(format(as.Date(paste0("2020-", start_date)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2020-", end_date)), "%j"))
  window_label <- paste0("Window (DOY ", start_doy, "-", end_doy, ")")
  
  # Convert lists to data frames if needed
  if (is.list(window_result) && !is.data.frame(window_result)) {
    window_result <- data.frame(
      mgmt_river = window_result$mgmt_river,
      total_production = window_result$total_production,
      percent_of_total_run = window_result$percent_of_total_run
    )
  }
  
  if (is.list(total_result) && !is.data.frame(total_result)) {
    total_result <- data.frame(
      mgmt_river = total_result$mgmt_river,
      total_production = total_result$total_production,
      percent_of_total_run = total_result$percent_of_total_run
    )
  }
  
  # Add percentages to edges for mapping
  edges$percent_of_total_run <- NA
  for (i in 1:nrow(window_result)) {
    mgmt_name <- window_result$mgmt_river[i]
    pct <- window_result$percent_of_total_run[i]
    edges$percent_of_total_run[edges$mgmt_river == mgmt_name] <- pct
  }
  
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  # Create map components
  map_filename <- paste0(watershed, "_", year, "_Window_", start_doy, "_", end_doy, "_Management.png")
  gg_hist <- create_doy_histogram(natal_data, window_data, window_label)
  
  # Main map
  main_plot <- ggplot() +
    geom_sf(data = basin, fill = "gray90", color = "gray70", size = 0.5) +
    geom_sf(data = managed_edges, aes(color = percent_of_total_run), size = 1.2) +
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "% Total Run",
                          na.value = "grey60", labels = function(x) paste0(round(x, 1), "%")) +
    coord_sf(datum = NA) +
    labs(title = paste0(window_label, " - ", watershed, " Management Units"),
         subtitle = paste("Year", year, "- Protection:", round(overall_pct, 1), "%")) +
    theme_void() + theme(legend.position = "right")
  
  # Bar chart - recreate the original comparison
  mgmt_comparison <- total_result %>%
    left_join(window_result %>% select(mgmt_river, percent_of_total_run) %>%
                rename(window_percent = percent_of_total_run), by = "mgmt_river") %>%
    mutate(window_percent = ifelse(is.na(window_percent), 0, window_percent)) %>%
    arrange(desc(percent_of_total_run))
  
  bar_plot <- ggplot(mgmt_comparison, aes(x = reorder(mgmt_river, percent_of_total_run))) +
    geom_col(aes(y = percent_of_total_run), fill = "gray80", alpha = 0.9) +
    geom_col(aes(y = window_percent), fill = "#8EB897", alpha = 0.9) +
    coord_flip() + scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
    labs(title = "Management Unit Contribution", x = "", y = "% Total Run") +
    theme_minimal() + theme(axis.text.y = element_text(size = 8))
  
  # Save combined plot using the original layout
  png(file.path(map_output_dir, map_filename), width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Create the layout properly without double-rendering
  if (!is.null(gg_hist)) {
    # If histogram exists, create 3-panel layout
    top_row <- arrangeGrob(main_plot, bar_plot, ncol = 2, widths = c(0.6, 0.4))
    final_plot <- grid.arrange(
      top_row, gg_hist,
      nrow = 2, heights = c(0.7, 0.3),
      top = paste0(window_label, " - ", watershed, " Management Units")
    )
  } else {
    # If no histogram, just create 2-panel layout
    final_plot <- grid.arrange(
      main_plot, bar_plot,
      ncol = 2, widths = c(0.6, 0.4),
      top = paste0(window_label, " - ", watershed, " Management Units")
    )
  }
  
  dev.off()
  message(paste("Created map:", map_filename))
}

# Helper function to create DOY histogram (from original code)
create_doy_histogram <- function(natal_data, window_data, window_label) {
  if (is.null(natal_data) || nrow(natal_data) == 0) return(NULL)
  
  tryCatch({
    ggplot() +
      geom_histogram(data = natal_data, aes(x = DOY), bins = 50, 
                     fill = "gray80", alpha = 0.7) +
      geom_histogram(data = window_data, aes(x = DOY), bins = 50, 
                     fill = "#8EB897", alpha = 0.8) +
      labs(title = paste("Fish Timing Distribution -", window_label),
           x = "Day of Year", y = "Count") +
      theme_minimal()
  }, error = function(e) {
    message("Could not create histogram: ", e$message)
    return(NULL)
  })
}

################################################################################
#                        SIMPLE FRONT-END CLOSURE ANALYSIS                    #
################################################################################
# PURPOSE: Calculate what % of salmon run is protected during closure windows
# INPUT:   Years, watershed, closure dates
# OUTPUT:  CSV with protection percentages by management unit
################################################################################

library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(here)
library(grid)
library(gridExtra)

# Source your existing functions (keep these as-is)
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

################################################################################
#                            MAIN ANALYSIS FUNCTION                           #
################################################################################

analyze_closure_simple <- function(years, 
                                   watershed = "Kusko",
                                   closure_start = "06-01", 
                                   closure_end = "06-11",
                                   create_maps = TRUE,
                                   clear_maps = TRUE) {
  
  # Convert dates to day-of-year
  start_doy <- as.numeric(format(as.Date(paste0("2024-", closure_start)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2024-", closure_end)), "%j"))
  
  message(paste("Analyzing closure window:", closure_start, "to", closure_end, 
                "(DOY", start_doy, "to", end_doy, ")"))
  
  # Set watershed parameters
  if (watershed == "Kusko") {
    params <- list(sensitivity_threshold = 0.7, min_error = 0.0006, min_stream_order = 3)
  } else if (watershed == "Yukon") {
    params <- list(sensitivity_threshold = 0.7, min_error = 0.003, min_stream_order = 5)
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon'")
  }
  
  all_results <- data.frame()
  
  # Setup map output directory if needed
  if (create_maps) {
    map_dir <- here(paste0("Basin Maps/FrontEnd_Closure_Management"))
    dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Clear existing maps if requested
    if (clear_maps) {
      existing_maps <- list.files(map_dir, pattern = "\\.png$", full.names = TRUE)
      if (length(existing_maps) > 0) {
        file.remove(existing_maps)
        message(paste("Cleared", length(existing_maps), "existing maps from:", map_dir))
      }
    }
    
    message(paste("Maps will be saved to:", map_dir))
  }
  
  ################################################################################
  #                              PROCESS EACH YEAR                              #
  ################################################################################
  
  for (year in years) {
    message(paste("\n--- Processing", watershed, "watershed for year", year, "---"))
    
    #--------------------------------------------------------------------------#
    #                           1. LOAD ALL DATA                              #
    #--------------------------------------------------------------------------#
    
    # Load spatial data (streams, basins, management units)
    spatial_data <- load_spatial_data(watershed, 8, params$min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    
    # Load fish data for this year
    natal_data <- load_natal_data(year, watershed)
    message(paste("Loaded", nrow(natal_data), "fish observations"))
    
    # Check for management units
    if (!"mgmt_river" %in% colnames(edges)) {
      stop("Management river data not found in edges. Use joined shapefile.")
    }
    
    #--------------------------------------------------------------------------#
    #                      2. FILTER TO CLOSURE WINDOW                        #
    #--------------------------------------------------------------------------#
    
    window_data <- natal_data %>% 
      filter(DOY >= start_doy & DOY <= end_doy)
    
    if (nrow(window_data) == 0) {
      warning(paste("No fish data in closure window for", year))
      next
    }
    
    message(paste("Found", nrow(window_data), "fish in closure window"))
    
    #--------------------------------------------------------------------------#
    #                       3. SETUP ASSIGNMENT PARAMETERS                    #
    #--------------------------------------------------------------------------#
    
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, params$min_error)
    priors <- setup_watershed_priors(edges, params$min_stream_order, watershed, natal_data)
    
    #--------------------------------------------------------------------------#
    #                    4. ASSIGN FISH TO NATAL STREAMS                      #
    #--------------------------------------------------------------------------#
    
    # Assign ALL fish to streams (full year)
    message("Assigning all fish to natal streams...")
    total_assignment <- perform_assignment(
      natal_data, edges, watershed, priors, pid_iso, error, params$sensitivity_threshold
    )
    total_stream_production <- apply(total_assignment, 1, sum, na.rm = TRUE)
    
    # Assign WINDOW fish to streams (closure period only)
    message("Assigning window fish to natal streams...")
    window_assignment <- perform_assignment(
      window_data, edges, watershed, priors, pid_iso, error, params$sensitivity_threshold
    )
    window_stream_production <- apply(window_assignment, 1, sum, na.rm = TRUE)
    
    #--------------------------------------------------------------------------#
    #                   5. AGGREGATE BY MANAGEMENT UNITS                      #
    #--------------------------------------------------------------------------#
    
    # Get total production by management unit (full year)
    total_mgmt <- process_mgmt_river_data(edges, total_stream_production)
    if (is.null(total_mgmt)) {
      warning(paste("No management unit data for", year))
      next
    }
    
    # Get window production by management unit (closure period)
    window_mgmt <- process_mgmt_river_data(edges, window_stream_production)
    
    # Remove spatial geometry if present
    if (inherits(total_mgmt, "sf")) total_mgmt <- st_drop_geometry(total_mgmt)
    if (inherits(window_mgmt, "sf")) window_mgmt <- st_drop_geometry(window_mgmt)
    
    #--------------------------------------------------------------------------#
    #                        6. CALCULATE PERCENTAGES                         #
    #--------------------------------------------------------------------------#
    
    # Total production across ALL management units for this year
    total_production_all_units <- sum(total_mgmt$total_production, na.rm = TRUE)
    
    # Join total and window data
    year_results <- total_mgmt %>%
      select(mgmt_river, total_production) %>%
      rename(total_annual_production = total_production) %>%
      left_join(
        window_mgmt %>% select(mgmt_river, total_production) %>%
          rename(window_production = total_production),
        by = "mgmt_river"
      ) %>%
      mutate(
        # Fill missing window production with 0
        window_production = ifelse(is.na(window_production), 0, window_production),
        
        # METRIC 1: What % of total run comes from this unit? (should sum to 100%)
        unit_contribution_pct = (total_annual_production / total_production_all_units) * 100,
        
        # METRIC 2: What % of this unit's fish are in the closure window?
        unit_protection_pct = ifelse(total_annual_production > 0,
                                     (window_production / total_annual_production) * 100, 0),
        
        # Add year info
        year = year
      ) %>%
      arrange(desc(unit_contribution_pct))
    
    all_results <- bind_rows(all_results, year_results)
    
    #--------------------------------------------------------------------------#
    #                          7. SUMMARY STATS                               #
    #--------------------------------------------------------------------------#
    
    total_window_protection <- sum(year_results$window_production, na.rm = TRUE)
    overall_protection_pct <- (total_window_protection / total_production_all_units) * 100
    
    message(paste("Overall protection for", year, ":", round(overall_protection_pct, 1), "%"))
    message(paste("Top 3 contributing units:"))
    print(year_results %>% slice_head(n = 3) %>% 
            select(mgmt_river, unit_contribution_pct, unit_protection_pct))
    
    #--------------------------------------------------------------------------#
    #                             8. CREATE MAPS                              #
    #--------------------------------------------------------------------------#
    
    if (create_maps) {
      # Use the EXACT original map function
      create_mgmt_closure_map(
        window_result = list(
          mgmt_river = year_results$mgmt_river,
          total_production = year_results$window_production,
          percent_of_total_run = (year_results$window_production / total_production_all_units) * 100
        ),
        total_result = list(
          mgmt_river = year_results$mgmt_river,
          total_production = year_results$total_annual_production,
          percent_of_total_run = year_results$unit_contribution_pct
        ),
        edges = edges,
        basin = basin,
        year = year,
        watershed = watershed,
        start_date = closure_start,
        end_date = closure_end,
        overall_pct = overall_protection_pct,
        map_output_dir = map_dir,
        natal_data = natal_data,
        window_data = window_data,
        total_production_all = total_production_all_units
      )
    }
  }
  
  ################################################################################
  #                               SAVE RESULTS                                  #
  ################################################################################
  
  if (nrow(all_results) == 0) {
    warning("No results to save")
    return(NULL)
  }
  
  # Create clean output
  final_results <- all_results %>%
    select(year, mgmt_river, unit_contribution_pct, unit_protection_pct) %>%
    rename(
      Year = year,
      Management_Unit = mgmt_river,
      Total_Contribution_Percent = unit_contribution_pct,
      Percent_Of_Unit_In_Window = unit_protection_pct
    )
  
  # Save to CSV
  output_dir <- here("Data/Frontend_Closure_Simple")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  filename <- paste0("closure_summary_", watershed, "_", 
                     min(years), "_to_", max(years), "_", 
                     closure_start, "_to_", closure_end, ".csv")
  
  write.csv(final_results, file.path(output_dir, filename), row.names = FALSE)
  message(paste("\n✓ Results saved to:", file.path(output_dir, filename)))
  
  ################################################################################
  #                              VERIFICATION                                    #
  ################################################################################
  
  # Check that unit contributions sum to ~100% each year
  verification <- final_results %>%
    group_by(Year) %>%
    summarise(
      total_contribution = sum(Total_Contribution_Percent, na.rm = TRUE),
      n_units = n(),
      .groups = "drop"
    )
  
  message("\n--- VERIFICATION ---")
  message("Unit contributions should sum to ~100% each year:")
  print(verification)
  
  if (all(abs(verification$total_contribution - 100) < 1)) {
    message("✓ All years sum correctly!")
  } else {
    warning("⚠ Some years don't sum to 100% - check calculations!")
  }
  
  return(final_results)
}

################################################################################
#                              EXAMPLE USAGE                                  #
################################################################################

# Run the analysis with maps (same as original)
results <- analyze_closure_simple(
  years = c(2017, 2018, 2019, 2020, 2021),
  watershed = "Kusko",
  closure_start = "06-01",
  closure_end = "06-11",
  create_maps = TRUE,   # Set to FALSE if you don't want maps
  clear_maps = TRUE     # Set to FALSE to keep existing maps
)

# View results
if (!is.null(results)) {
  message("\n--- SAMPLE RESULTS ---")
  print(head(results, 10))
}