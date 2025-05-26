# simple_closure_analysis.R
# Streamlined front-end closure window analysis using EXACT original calculation method

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

# ==============================================================================
# MAIN CLOSURE ANALYSIS FUNCTION (using exact original calculation method)
# ==============================================================================

#' Analyze front-end closure window using exact original method from Calculate_Huc_Closure_protection.R
#'
#' @param years Vector of years to analyze
#' @param watershed "Kusko" or "Yukon" 
#' @param closure_start Start date "MM-DD" format (default "06-01")
#' @param closure_end End date "MM-DD" format (default "06-11")
#' @param create_maps Whether to create basin maps (default TRUE)
#' @return Data frame with results
analyze_closure_window <- function(years, 
                                   watershed = "Kusko",
                                   closure_start = "06-01", 
                                   closure_end = "06-11",
                                   create_maps = TRUE) {
  
  # Create output directories
  csv_dir <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Frontend_Closure"
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (create_maps) {
    map_output_dir <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Basin Maps/FrontEnd_Closure"
    dir.create(map_output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Set watershed parameters (exactly like original)
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
  
  # Convert closure dates to DOY (exactly like original)
  closure_start_doy <- as.numeric(format(as.Date(paste0("2024-", closure_start)), "%j"))
  closure_end_doy <- as.numeric(format(as.Date(paste0("2024-", closure_end)), "%j"))
  
  message(paste("Analyzing closure window DOY", closure_start_doy, "to", closure_end_doy, 
                "(", closure_start, "to", closure_end, ")"))
  
  # Initialize results
  all_results <- data.frame()
  
  # Process each year (exactly like original Calculate_Huc_Closure_protection.R)
  for (year in years) {
    message(paste("Processing", watershed, "watershed for year", year))
    
    # Load spatial data (exactly like original)
    spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    Huc <- spatial_data$Huc
    
    # Load natal origins data (exactly like original)
    natal_data <- load_natal_data(year, watershed)
    
    # Filter the natal data to only include fish within the closure window (exactly like original)
    window_data <- natal_data %>% 
      filter(DOY >= closure_start_doy & DOY <= closure_end_doy)
    
    if (nrow(window_data) == 0) {
      warning(paste("No data found in closure window for", watershed, "in", year))
      next
    }
    
    message(paste("  Found", nrow(window_data), "of", nrow(natal_data), 
                  "total observations within closure window"))
    
    # Extract isoscape prediction and error values (exactly like original)
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    
    # Calculate error values (exactly like original)
    error <- calculate_error(pid_isose, min_error)
    
    # Set up watershed-specific priors (exactly like original)
    priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
    
    # --- Calculate total production across all DOY values (exactly like original) ---
    message("  Calculating total production across all DOY values...")
    
    total_assignment_matrix <- perform_assignment(
      natal_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Process total assignments to get basin-scale values (exactly like original)
    total_basin_results <- process_assignments(total_assignment_matrix)
    total_basin_assign_sum <- total_basin_results$sum
    
    # Process HUC data for total production (exactly like original)
    total_huc_result <- process_huc_data(edges, basin, Huc, total_basin_assign_sum, 8)
    
    # --- Calculate production within closure window (exactly like original) ---
    message("  Calculating production within closure window...")
    
    window_assignment_matrix <- perform_assignment(
      window_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Process window assignments to get basin-scale values (exactly like original)
    window_basin_results <- process_assignments(window_assignment_matrix)
    window_basin_assign_sum <- window_basin_results$sum
    
    # Process HUC data for window production (exactly like original)
    window_huc_result <- process_huc_data(edges, basin, Huc, window_basin_assign_sum, 8)
    
    # Extract necessary data from both results (exactly like original)
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
    
    # Calculate total production for normalization (exactly like original)
    total_production_all_hucs <- sum(total_huc_df$total_production, na.rm = TRUE)
    
    # Calculate the overall protection percentage for the entire run (exactly like original)
    overall_protection_pct <- sum(window_huc_df$total_production, na.rm = TRUE) / 
      max(total_production_all_hucs, 0.0001) * 100  # Avoid division by zero
    
    # Calculate CPUE proportions (like original)
    total_cpue <- sum(natal_data$dailyCPUEprop, na.rm = TRUE)
    window_cpue <- sum(window_data$dailyCPUEprop, na.rm = TRUE)
    window_cpue_proportion <- window_cpue / total_cpue
    
    # Join the datasets (exactly like original)
    huc_col <- "HUC8"
    name_col <- "Name"
    
    # Use safer join approach (exactly like original)
    combined_results <- total_huc_df %>%
      select(all_of(c(huc_col, name_col, "total_production"))) %>%
      rename(total_production_all = total_production) 
    
    # Add window production (exactly like original)
    window_production <- window_huc_df %>%
      select(all_of(c(huc_col, "total_production"))) %>%
      rename(total_production_window = total_production)
    
    combined_results <- combined_results %>%
      left_join(window_production, by = huc_col)
    
    # Replace NA values with 0 (exactly like original)
    combined_results$total_production_window[is.na(combined_results$total_production_window)] <- 0
    
    # Calculate percentages (exactly like original)
    combined_results <- combined_results %>%
      mutate(
        # Calculate total percentages (gray bars + green bars)
        huc_pct_of_total_run = (total_production_all / max(total_production_all_hucs, 0.0001)) * 100,
        
        # Calculate window percentages (green bars only)
        huc_pct_in_window = (total_production_window / max(total_production_all_hucs, 0.0001)) * 100,
        
        # Calculate the protection percentage (green as % of that HUC's total)
        protection_percentage = ifelse(total_production_all > 0,
                                       (total_production_window / total_production_all) * 100,
                                       0)
      ) %>%
      # Add year and metadata information (exactly like original)
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
      # Sort by HUC percentage of total run in descending order (exactly like original)
      arrange(desc(huc_pct_of_total_run))
    
    # Add to overall results
    all_results <- bind_rows(all_results, combined_results)
    
    # Create basin map if requested (using exact original method)
    if (create_maps) {
      # Calculate the correct percent_of_total_run values for the map
      # This should be the window production as % of total annual production
      window_huc_result_for_map <- window_huc_result
      window_huc_result_for_map$percent_of_total_run <- (window_huc_result_for_map$total_production / total_production_all_hucs) * 100
      
      # Also add this to total_huc_result for the comparison
      total_huc_result_for_map <- total_huc_result
      total_huc_result_for_map$percent_of_total_run <- (total_huc_result_for_map$total_production / total_production_all_hucs) * 100
      
      create_closure_map(window_huc_result_for_map, total_huc_result_for_map, year, watershed, 
                         closure_start, closure_end, overall_protection_pct, 
                         map_output_dir, natal_data, window_data)
    }
  }
  
  # Check if we have results
  if (nrow(all_results) == 0) {
    warning("No results collected - unable to save CSV")
    return(NULL)
  }
  
  # Export one simplified CSV with all years combined
  csv_filename <- paste0("huc_closure_summary_", watershed, "_", 
                         min(years), "_to_", max(years), "_", 
                         closure_start, "_to_", closure_end, ".csv")
  csv_path <- file.path(csv_dir, csv_filename)
  
  # Create simplified results with just the 3 key columns for all years
  simplified_results <- all_results %>%
    select(year, Name, huc_pct_of_total_run, protection_percentage) %>%
    rename(
      "Year" = year,
      "HUC_Name" = Name,
      "Total_HUC_Contribution_Percent" = huc_pct_of_total_run,
      "Proportion_in_Closure_Window_Percent" = protection_percentage
    )
  
  write.csv(simplified_results, csv_path, row.names = FALSE)
  message(paste("Saved simplified results with all years to:", csv_path))
  
  return(simplified_results)
}

# ==============================================================================
# MAP CREATION FUNCTION (using exact original method from fishing_window_analysis.R)
# ==============================================================================

#' Create basin map exactly like original fishing_window_analysis.R
create_closure_map <- function(window_huc_result, all_hucs_result, year, watershed, 
                               start_date, end_date, overall_pct, map_output_dir, 
                               natal_data, window_data) {
  
  # Convert dates to DOY for title
  start_doy <- as.numeric(format(as.Date(paste0("2020-", start_date)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2020-", end_date)), "%j"))
  
  # Create window label exactly like original
  window_label <- paste0("Fishing Window (DOY ", start_doy, "-", end_doy, ", ", 
                         format(as.Date(paste0("2020-", start_date)), "%b %d"), "-",
                         format(as.Date(paste0("2020-", end_date)), "%b %d"), ")")
  
  # Use exact same filename as original
  map_filename <- paste0(watershed, "_", year, "_FishingWindow_", start_doy, "_", end_doy, "_HUC8.png")
  map_path <- file.path(map_output_dir, map_filename)
  
  # Create improved histogram highlighting the fishing window (exactly like original)
  gg_hist <- create_doy_histogram(natal_data, window_data, window_label)
  
  # Create PNG with exact original specifications
  png(file = map_path, width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Set up the plotting layout (exact same as original)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, 
                                                               heights = grid::unit(c(0.7, 0.3), "npc"),
                                                               widths = grid::unit(c(0.6, 0.4), "npc"))))
  
  # Create main map plot with percent of total run (exact original styling)
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
      subtitle = paste("Year", year, "- Overall Window Percentage:", 
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
  
  # Create the histogram showing total vs. window contribution (exactly like original)
  huc_histogram <- create_huc_total_quartile_histogram(
    final_result = window_huc_result,
    all_hucs_data = all_hucs_result,
    quartile_label = window_label
  )
  
  # Print the plots to the PNG (exact same layout as original)
  print(main_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(huc_histogram, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  # Add the DOY histogram at the bottom spanning both columns (exactly like original)
  if (!is.null(gg_hist)) {
    gg_hist <- enforce_histogram_limits(gg_hist) + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      # Add vertical lines at the window boundaries (exactly like original)
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
  
  message(paste("Map saved:", map_filename))
}

# ==============================================================================
# USAGE
# ==============================================================================

# Run analysis for Kuskokwim 2017-2021 (exactly like original)
results <- analyze_closure_window(
  years = c(2017, 2018, 2019, 2020, 2021),
  watershed = "Kusko",
  closure_start = "06-01",
  closure_end = "06-11",
  create_maps = TRUE
)