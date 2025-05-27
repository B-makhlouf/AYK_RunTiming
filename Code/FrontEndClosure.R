################################################################################
# FRONTENDCLOSURE.R - FRONT-END CLOSURE WINDOW ANALYSIS
################################################################################
# PURPOSE: Analyzes what proportion of salmon run is protected during closure windows
# OUTPUT: Shows percentage of total run protected during specified closure dates
# TYPICAL USE: June 1-11 closure window for Kuskokwim management
################################################################################

library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(here)
library(grid)
library(gridExtra)

# Source existing utility functions
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

################################################################################
# MAIN CLOSURE ANALYSIS FUNCTION
################################################################################

#' Analyze front-end closure window protection effectiveness
#' 
#' @param years Vector of years to analyze
#' @param watershed "Kusko" or "Yukon" 
#' @param closure_start Start date in "MM-DD" format (default "06-01")
#' @param closure_end End date in "MM-DD" format (default "06-11")
#' @param create_maps Whether to create basin maps (default TRUE)
#' @return Data frame with results showing protection percentages by HUC
analyze_closure_window <- function(years, 
                                   watershed = "Kusko",
                                   closure_start = "06-01", 
                                   closure_end = "06-11",
                                   create_maps = TRUE) {
  
  ################################################################################
  # SETUP: CREATE DIRECTORIES AND SET PARAMETERS
  ################################################################################
  
  # Create output directories
  csv_dir <- here("Data/Frontend_Closure")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (create_maps) {
    map_output_dir <- here("Basin Maps/FrontEnd_Closure")
    dir.create(map_output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Set watershed-specific parameters
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
  
  ################################################################################
  # CONVERT CLOSURE DATES TO DAY OF YEAR
  ################################################################################
  
  closure_start_doy <- as.numeric(format(as.Date(paste0("2024-", closure_start)), "%j"))
  closure_end_doy <- as.numeric(format(as.Date(paste0("2024-", closure_end)), "%j"))
  
  message(paste("Analyzing closure window DOY", closure_start_doy, "to", closure_end_doy, 
                "(", closure_start, "to", closure_end, ")"))
  
  # Initialize results storage
  all_results <- data.frame()
  
  ################################################################################
  # PROCESS EACH YEAR
  ################################################################################
  
  for (year in years) {
    message(paste("Processing", watershed, "watershed for year", year))
    
    ################################################################################
    # LOAD DATA FOR THIS YEAR
    ################################################################################
    
    # Load spatial data
    spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    Huc <- spatial_data$Huc
    
    # Load natal origins data
    natal_data <- load_natal_data(year, watershed)
    
    ################################################################################
    # FILTER DATA TO CLOSURE WINDOW
    ################################################################################
    
    # Filter to only fish caught within the closure window
    window_data <- natal_data %>% 
      filter(DOY >= closure_start_doy & DOY <= closure_end_doy)
    
    if (nrow(window_data) == 0) {
      warning(paste("No data found in closure window for", watershed, "in", year))
      next
    }
    
    message(paste("Found", nrow(window_data), "of", nrow(natal_data), 
                  "total observations within closure window"))
    
    ################################################################################
    # SETUP ASSIGNMENT PARAMETERS
    ################################################################################
    
    # Extract isoscape prediction and error values
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, min_error)
    priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
    
    ################################################################################
    # CALCULATE TOTAL ANNUAL PRODUCTION (FOR NORMALIZATION)
    ################################################################################
    
    message("Calculating total production across all DOY values...")
    
    total_assignment_matrix <- perform_assignment(
      natal_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    total_basin_results <- process_assignments(total_assignment_matrix)
    total_basin_assign_sum <- total_basin_results$sum
    total_huc_result <- process_huc_data(edges, basin, Huc, total_basin_assign_sum, 8)
    
    ################################################################################
    # CALCULATE PRODUCTION WITHIN CLOSURE WINDOW
    ################################################################################
    
    message("Calculating production within closure window...")
    
    window_assignment_matrix <- perform_assignment(
      window_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    window_basin_results <- process_assignments(window_assignment_matrix)
    window_basin_assign_sum <- window_basin_results$sum
    window_huc_result <- process_huc_data(edges, basin, Huc, window_basin_assign_sum, 8)
    
    ################################################################################
    # CALCULATE PROTECTION PERCENTAGES
    ################################################################################
    
    # Extract data from results
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
    
    # Calculate overall protection percentage for entire run
    overall_protection_pct <- sum(window_huc_df$total_production, na.rm = TRUE) / 
      max(total_production_all_hucs, 0.0001) * 100
    
    # Calculate CPUE proportions
    total_cpue <- sum(natal_data$dailyCPUEprop, na.rm = TRUE)
    window_cpue <- sum(window_data$dailyCPUEprop, na.rm = TRUE)
    window_cpue_proportion <- window_cpue / total_cpue
    
    ################################################################################
    # COMBINE TOTAL AND WINDOW RESULTS BY HUC
    ################################################################################
    
    huc_col <- "HUC8"
    name_col <- "Name"
    
    # Combine total and window production data
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
    
    ################################################################################
    # CALCULATE FINAL PROTECTION METRICS
    ################################################################################
    
    combined_results <- combined_results %>%
      mutate(
        # Total HUC contribution to entire run
        huc_pct_of_total_run = (total_production_all / max(total_production_all_hucs, 0.0001)) * 100,
        
        # HUC contribution during closure window
        huc_pct_in_window = (total_production_window / max(total_production_all_hucs, 0.0001)) * 100,
        
        # Protection percentage (window production as % of that HUC's total)
        protection_percentage = ifelse(total_production_all > 0,
                                       (total_production_window / total_production_all) * 100,
                                       0)
      ) %>%
      # Add metadata
      mutate(
        year = year,
        watershed = watershed,
        closure_start_date = closure_start,
        closure_end_date = closure_end,
        closure_start_doy = closure_start_doy,
        closure_end_doy = closure_end_doy,
        overall_window_production_pct = overall_protection_pct,
        overall_window_cpue_pct = window_cpue_proportion * 100
      ) %>%
      # Sort by HUC percentage of total run (most important HUCs first)
      arrange(desc(huc_pct_of_total_run))
    
    # Add to overall results
    all_results <- bind_rows(all_results, combined_results)
    
    ################################################################################
    # CREATE BASIN MAP IF REQUESTED
    ################################################################################
    
    if (create_maps) {
      # Prepare HUC results for mapping
      window_huc_result_for_map <- window_huc_result
      window_huc_result_for_map$percent_of_total_run <- (window_huc_result_for_map$total_production / total_production_all_hucs) * 100
      
      total_huc_result_for_map <- total_huc_result
      total_huc_result_for_map$percent_of_total_run <- (total_huc_result_for_map$total_production / total_production_all_hucs) * 100
      
      create_closure_map(window_huc_result_for_map, total_huc_result_for_map, year, watershed, 
                         closure_start, closure_end, overall_protection_pct, 
                         map_output_dir, natal_data, window_data)
    }
  }
  
  ################################################################################
  # EXPORT RESULTS TO CSV
  ################################################################################
  
  if (nrow(all_results) == 0) {
    warning("No results collected - unable to save CSV")
    return(NULL)
  }
  
  # Create simplified CSV with key results
  csv_filename <- paste0("huc_closure_summary_", watershed, "_", 
                         min(years), "_to_", max(years), "_", 
                         closure_start, "_to_", closure_end, ".csv")
  csv_path <- file.path(csv_dir, csv_filename)
  
  # Simplified results with key columns
  simplified_results <- all_results %>%
    select(year, Name, huc_pct_of_total_run, protection_percentage) %>%
    rename(
      "Year" = year,
      "HUC_Name" = Name,
      "Total_HUC_Contribution_Percent" = huc_pct_of_total_run,
      "Proportion_in_Closure_Window_Percent" = protection_percentage
    )
  
  write.csv(simplified_results, csv_path, row.names = FALSE)
  message(paste("Saved simplified results to:", csv_path))
  
  return(simplified_results)
}

################################################################################
# MAP CREATION FUNCTION
################################################################################

#' Create closure protection map showing window effectiveness
create_closure_map <- function(window_huc_result, all_hucs_result, year, watershed, 
                               start_date, end_date, overall_pct, map_output_dir, 
                               natal_data, window_data) {
  
  # Convert dates to DOY for title
  start_doy <- as.numeric(format(as.Date(paste0("2020-", start_date)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2020-", end_date)), "%j"))
  
  # Create descriptive window label
  window_label <- paste0("Fishing Window (DOY ", start_doy, "-", end_doy, ", ", 
                         format(as.Date(paste0("2020-", start_date)), "%b %d"), "-",
                         format(as.Date(paste0("2020-", end_date)), "%b %d"), ")")
  
  # Create filename
  map_filename <- paste0(watershed, "_", year, "_FishingWindow_", start_doy, "_", end_doy, "_HUC8.png")
  map_path <- file.path(map_output_dir, map_filename)
  
  # Create histogram highlighting the closure window
  gg_hist <- create_doy_histogram(natal_data, window_data, window_label)
  
  ################################################################################
  # CREATE PNG MAP
  ################################################################################
  
  png(file = map_path, width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Set up plotting layout
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, 
                                                               heights = grid::unit(c(0.7, 0.3), "npc"),
                                                               widths = grid::unit(c(0.6, 0.4), "npc"))))
  
  # Main map plot showing percent of total run protected
  main_plot <- ggplot() +
    geom_sf(data = window_huc_result, aes(fill = percent_of_total_run), color = "white", size = 0.1) +
    scale_fill_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Percent of\nTotal Run",
      na.value = "grey95",
      labels = function(x) paste0(round(x, 1), "%"),
      guide = guide_colorbar(
        barwidth = 1, barheight = 15,
        frame.colour = "grey40", ticks.colour = "grey40",
        show.limits = TRUE
      )
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0(window_label, ": Percent of Total Run - ", watershed, " Watershed"),
      subtitle = paste("Year", year, "- Overall Window Protection:", 
                       round(overall_pct, 1), "%")
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "grey30"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold", color = "grey30"),
      legend.text = element_text(color = "grey30"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  # Comparison histogram showing total vs. window contribution
  huc_histogram <- create_huc_total_quartile_histogram(
    final_result = window_huc_result,
    all_hucs_data = all_hucs_result,
    quartile_label = window_label
  )
  
  # Plot main map and histogram
  print(main_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(huc_histogram, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  # Add DOY histogram at bottom with closure window marked
  if (!is.null(gg_hist)) {
    gg_hist <- enforce_histogram_limits(gg_hist) + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      # Add vertical lines marking closure window boundaries
      geom_vline(xintercept = start_doy, color = "blue", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = end_doy, color = "blue", linetype = "dashed", linewidth = 1) +
      annotate("text", x = start_doy - 2, y = 0.09, 
               label = paste0("Start: DOY ", start_doy, "\n", format(as.Date(paste0("2020-", start_date)), "%b %d")),
               color = "blue", hjust = 1) +
      annotate("text", x = end_doy + 2, y = 0.09, 
               label = paste0("End: DOY ", end_doy, "\n", format(as.Date(paste0("2020-", end_date)), "%b %d")),
               color = "blue", hjust = 0)
    
    print(gg_hist, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  }
  
  dev.off()
  
  message(paste("Closure map saved:", map_filename))
}

################################################################################
# EXAMPLE USAGE (COMMENTED OUT)
################################################################################

# Run closure analysis for Kuskokwim 2017-2021 with June 1-11 closure
# results <- analyze_closure_window(
#   years = c(2017, 2018, 2019, 2020, 2021),
#   watershed = "Kusko",
#   closure_start = "06-01",
#   closure_end = "06-11",
#   create_maps = TRUE
# )