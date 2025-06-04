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
#' @param use_mgmt_units Whether to use management units instead of HUCs (default TRUE)
#' @return Data frame with results showing protection percentages by unit
analyze_closure_window <- function(years, 
                                   watershed = "Kusko",
                                   closure_start = "06-01", 
                                   closure_end = "06-11",
                                   create_maps = TRUE,
                                   use_mgmt_units = TRUE) {
  
  # Setup directories and parameters
  unit_type <- ifelse(use_mgmt_units, "Management", "HUC")
  csv_dir <- here(paste0("Data/Frontend_Closure_", unit_type))
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (create_maps) {
    map_output_dir <- here(paste0("Basin Maps/FrontEnd_Closure_", unit_type))
    dir.create(map_output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Watershed parameters
  params <- if (watershed == "Yukon") {
    list(sensitivity_threshold = 0.7, min_error = 0.003, min_stream_order = 5)
  } else if (watershed == "Kusko") {
    list(sensitivity_threshold = 0.7, min_error = 0.0006, min_stream_order = 3)
  } else {
    stop(paste("Unknown watershed:", watershed))
  }
  
  # Convert closure dates to DOY
  closure_start_doy <- as.numeric(format(as.Date(paste0("2024-", closure_start)), "%j"))
  closure_end_doy <- as.numeric(format(as.Date(paste0("2024-", closure_end)), "%j"))
  
  message(paste("Analyzing closure window DOY", closure_start_doy, "to", closure_end_doy, 
                "using", unit_type, "units"))
  
  all_results <- data.frame()
  
  ################################################################################
  # PROCESS EACH YEAR
  ################################################################################
  
  for (year in years) {
    message(paste("Processing", watershed, "watershed for year", year))
    
    # Load data
    spatial_data <- load_spatial_data(watershed, 8, params$min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    Huc <- spatial_data$Huc
    natal_data <- load_natal_data(year, watershed)
    
    # Check for management units if requested
    if (use_mgmt_units && !"mgmt_river" %in% colnames(edges)) {
      stop("Management river data not found. Set use_mgmt_units = FALSE or use joined shapefile.")
    }
    
    # Filter to closure window
    window_data <- natal_data %>% 
      filter(DOY >= closure_start_doy & DOY <= closure_end_doy)
    
    if (nrow(window_data) == 0) {
      warning(paste("No data found in closure window for", watershed, "in", year))
      next
    }
    
    message(paste("Found", nrow(window_data), "of", nrow(natal_data), 
                  "observations within closure window"))
    
    # Setup assignment parameters
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, params$min_error)
    priors <- setup_watershed_priors(edges, params$min_stream_order, watershed, natal_data)
    
    # Calculate total annual production
    total_assignment_matrix <- perform_assignment(
      natal_data, edges, watershed, priors, pid_iso, error, params$sensitivity_threshold)
    total_basin_assign_sum <- apply(total_assignment_matrix, 1, sum, na.rm = TRUE)
    
    # Calculate window production
    window_assignment_matrix <- perform_assignment(
      window_data, edges, watershed, priors, pid_iso, error, params$sensitivity_threshold)
    window_basin_assign_sum <- apply(window_assignment_matrix, 1, sum, na.rm = TRUE)
    
    # Process by unit type
    if (use_mgmt_units) {
      total_result <- process_mgmt_river_data(edges, total_basin_assign_sum)
      window_result <- process_mgmt_river_data(edges, window_basin_assign_sum)
      
      if (is.null(total_result) || is.null(window_result)) {
        warning(paste("No management unit data for", watershed, "in", year))
        next
      }
      
      unit_col <- "mgmt_river"
      name_col <- "mgmt_river"
    } else {
      total_result <- process_huc_data(edges, basin, Huc, total_basin_assign_sum, 8)
      window_result <- process_huc_data(edges, basin, Huc, window_basin_assign_sum, 8)
      unit_col <- "HUC8"
      name_col <- "Name"
    }
    
    # Calculate metrics
    if (inherits(total_result, "sf")) total_result <- st_drop_geometry(total_result)
    if (inherits(window_result, "sf")) window_result <- st_drop_geometry(window_result)
    
    total_production_all <- sum(total_result$total_production, na.rm = TRUE)
    overall_protection_pct <- sum(window_result$total_production, na.rm = TRUE) / 
      max(total_production_all, 0.0001) * 100
    
    # Combine results
    combined_results <- total_result %>%
      select(all_of(c(name_col, "total_production"))) %>%
      rename(total_production_all = total_production) %>%
      left_join(
        window_result %>% select(all_of(c(name_col, "total_production"))) %>%
          rename(total_production_window = total_production),
        by = name_col
      ) %>%
      mutate(
        total_production_window = ifelse(is.na(total_production_window), 0, total_production_window),
        unit_pct_of_total_run = (total_production_all / max(total_production_all, 0.0001)) * 100,
        protection_percentage = ifelse(total_production_all > 0,
                                       (total_production_window / total_production_all) * 100, 0),
        year = year,
        watershed = watershed,
        closure_start_date = closure_start,
        closure_end_date = closure_end,
        overall_window_production_pct = overall_protection_pct
      ) %>%
      arrange(desc(unit_pct_of_total_run))
    
    all_results <- bind_rows(all_results, combined_results)
    
    # Create maps if requested
    if (create_maps) {
      if (use_mgmt_units) {
        create_mgmt_closure_map(window_result, total_result, edges, basin, year, watershed, 
                                closure_start, closure_end, overall_protection_pct, 
                                map_output_dir, natal_data, window_data, total_production_all)
      } else {
        create_closure_map(window_result, total_result, year, watershed, 
                           closure_start, closure_end, overall_protection_pct, 
                           map_output_dir, natal_data, window_data)
      }
    }
  }
  
  # Export results with both metrics
  if (nrow(all_results) == 0) {
    warning("No results collected")
    return(NULL)
  }
  
  csv_filename <- paste0(tolower(unit_type), "_closure_summary_", watershed, "_", 
                         min(years), "_to_", max(years), "_", 
                         closure_start, "_to_", closure_end, ".csv")
  
  simplified_results <- all_results %>%
    select(year, all_of(name_col), unit_pct_of_total_run, protection_percentage) %>%
    setNames(c("Year", "Management_Unit", 
               "Total_Contribution_Percent", 
               "Percent_Of_Unit_In_Window"))
  
  write.csv(simplified_results, file.path(csv_dir, csv_filename), row.names = FALSE)
  message(paste("Saved results to:", file.path(csv_dir, csv_filename)))
  
  return(simplified_results)
}

################################################################################
# MANAGEMENT UNIT MAP CREATION
################################################################################

create_mgmt_closure_map <- function(window_result, total_result, edges, basin, year, watershed, 
                                    start_date, end_date, overall_pct, map_output_dir, 
                                    natal_data, window_data, total_production_all) {
  
  start_doy <- as.numeric(format(as.Date(paste0("2020-", start_date)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2020-", end_date)), "%j"))
  window_label <- paste0("Window (DOY ", start_doy, "-", end_doy, ")")
  
  # Add percentages to results and edges
  window_result$percent_of_total_run <- (window_result$total_production / total_production_all) * 100
  total_result$percent_of_total_run <- (total_result$total_production / total_production_all) * 100
  
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
  
  # Bar chart
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
  
  # Save combined plot
  png(file.path(map_output_dir, map_filename), width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Use gridExtra for layout instead
  combined_plot <- grid.arrange(
    main_plot, bar_plot,
    ncol = 2, widths = c(0.6, 0.4),
    top = paste0(window_label, " - ", watershed, " Management Units")
  )
  
  if (!is.null(gg_hist)) {
    final_plot <- grid.arrange(
      combined_plot, gg_hist,
      nrow = 2, heights = c(0.7, 0.3)
    )
  }
  
  dev.off()
  message(paste("Created map:", map_filename))
}

################################################################################
# HUC MAP CREATION (ORIGINAL)
################################################################################

create_closure_map <- function(window_huc_result, all_hucs_result, year, watershed, 
                               start_date, end_date, overall_pct, map_output_dir, 
                               natal_data, window_data) {
  
  start_doy <- as.numeric(format(as.Date(paste0("2020-", start_date)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2020-", end_date)), "%j"))
  window_label <- paste0("Window (DOY ", start_doy, "-", end_doy, ")")
  
  map_filename <- paste0(watershed, "_", year, "_Window_", start_doy, "_", end_doy, "_HUC8.png")
  gg_hist <- create_doy_histogram(natal_data, window_data, window_label)
  
  main_plot <- ggplot() +
    geom_sf(data = window_huc_result, aes(fill = percent_of_total_run), color = "white", size = 0.1) +
    scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "% Total Run",
                         labels = function(x) paste0(round(x, 1), "%")) +
    coord_sf(datum = NA) +
    labs(title = paste0(window_label, " - ", watershed, " HUCs"),
         subtitle = paste("Year", year, "- Protection:", round(overall_pct, 1), "%")) +
    theme_void() + theme(legend.position = "right")
  
  huc_histogram <- create_huc_total_quartile_histogram(window_huc_result, all_hucs_result, window_label)
  
  png(file.path(map_output_dir, map_filename), width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  combined_plot <- grid.arrange(
    main_plot, huc_histogram,
    ncol = 2, widths = c(0.6, 0.4)
  )
  
  if (!is.null(gg_hist)) {
    final_plot <- grid.arrange(
      combined_plot, gg_hist,
      nrow = 2, heights = c(0.7, 0.3)
    )
  }
  
  dev.off()
  message(paste("Created map:", map_filename))
}

################################################################################
# EXAMPLE USAGE
################################################################################

# Run closure analysis with management units (default)
results_mgmt <- analyze_closure_window(
  years = c(2017, 2018, 2019, 2020, 2021),
  watershed = "Kusko",
  closure_start = "06-01",
  closure_end = "06-11",
  create_maps = TRUE,
  use_mgmt_units = TRUE
)

# Or run with HUCs
# results_huc <- analyze_closure_window(
#   years = c(2017, 2018, 2019, 2020, 2021),
#   watershed = "Kusko", 
#   closure_start = "06-01",
#   closure_end = "06-11",
#   create_maps = TRUE,
#   use_mgmt_units = FALSE
# )