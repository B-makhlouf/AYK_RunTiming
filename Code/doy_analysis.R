################################################################################
# DOY_ANALYSIS.R - DAY OF YEAR QUARTILE ANALYSIS
################################################################################
# PURPOSE: Analyzes salmon run timing by dividing data into 4 DOY quartiles
# OUTPUT: Shows proportion of production within each timing quartile
# KEY: Each quartile map shows what % of that quartile's production comes from each Management Unit
################################################################################

library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)

# Source required utilities
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

################################################################################
# MAIN DOY QUARTILE ANALYSIS FUNCTION
################################################################################

#' Perform DOY Quartile Analysis
#' Creates maps showing proportion of production within each timing quartile
DOY_Quartile_Analysis <- function(year, watershed, sensitivity_threshold, min_error, 
                                  min_stream_order = 3, 
                                  return_values = FALSE) {
  
  message(paste("=== Starting DOY Quartile Analysis for", year, watershed, "==="))
  
  ################################################################################
  # SETUP: LOAD DATA AND PARAMETERS
  ################################################################################
  
  # Load spatial and natal data
  spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  
  natal_data <- load_natal_data(year, watershed)
  message(paste("Loaded", nrow(natal_data), "natal data records"))
  
  # Divide data into DOY quartiles
  quartile_data <- divide_doy_quartiles(natal_data)
  quartile_subsets <- quartile_data$subsets
  subset_labels <- quartile_data$labels
  
  # Setup assignment parameters
  pid_iso <- edges$iso_pred
  pid_isose <- edges$isose_pred
  error <- calculate_error(pid_isose, min_error)
  priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
  
  ################################################################################
  # CREATE OUTPUT DIRECTORIES
  ################################################################################
  
  dir.create(here("Basin Maps/DOY_Quartile/Tribs"), showWarnings = FALSE, recursive = TRUE)
  
  # Storage for return values
  if (return_values) {
    all_results <- list()
  }
  
  ################################################################################
  # PROCESS EACH QUARTILE
  ################################################################################
  
  for (q in 1:length(quartile_subsets)) {
    current_subset <- quartile_subsets[[q]]
    
    # Skip empty quartiles
    if (nrow(current_subset) == 0) {
      message(paste("Skipping", subset_labels[q], "- no data"))
      next
    }
    
    message(paste("Processing", subset_labels[q], "with", nrow(current_subset), "data points"))
    
    ################################################################################
    # PERFORM ASSIGNMENT FOR THIS QUARTILE
    ################################################################################
    
    # Create unique ID for output files
    subset_id <- paste0(watershed, "_", year, "_DOY_Q", q)
    
    # Perform Bayesian assignment
    assignment_matrix <- perform_assignment(
      current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Process assignments to get basin-scale values
    basin_results <- process_assignments(assignment_matrix)
    basin_assign_rescale <- basin_results$rescale
    basin_assign_norm <- basin_results$norm
    
    ################################################################################
    # PROCESS MANAGEMENT RIVER DATA - SHOWS PROPORTION WITHIN THIS QUARTILE
    ################################################################################
    
    mgmt_result <- process_mgmt_river_data(edges, basin_assign_rescale)
    
    ################################################################################
    # CREATE VISUALIZATION PLOTS
    ################################################################################
    
    # Create DOY histogram showing this quartile highlighted
    gg_hist <- create_doy_histogram(natal_data, current_subset, subset_labels[q])
    
    ################################################################################
    # CREATE MANAGEMENT RIVER MAP - SHOWS PRODUCTION PROPORTION
    ################################################################################
    
    if (!is.null(mgmt_result)) {
      mgmt_filepath <- file.path(here("Basin Maps/DOY_Quartile/Management"), 
                                 paste0(subset_id, "_Management.png"))
      
      # Create output directory
      dir.create(here("Basin Maps/DOY_Quartile/Management"), showWarnings = FALSE, recursive = TRUE)
      
      # Create management river map
      create_simple_mgmt_map(mgmt_result, edges, basin, gg_hist, year, watershed, 
                             sensitivity_threshold, min_stream_order, subset_labels[q], mgmt_filepath)
    }
    
    ################################################################################
    # CREATE TRIBUTARY MAP
    ################################################################################
    
    trib_filepath <- file.path(here("Basin Maps/DOY_Quartile/Tribs"), 
                               paste0(subset_id, ".png"))
    
    create_doy_tributary_map(
      basin = basin,
      edges = edges,
      basin_assign_norm = basin_assign_norm,
      StreamOrderPrior = priors$StreamOrderPrior,
      pid_prior = priors$pid_prior,
      gg_hist = gg_hist,
      year = year,
      watershed = watershed,
      sensitivity_threshold = sensitivity_threshold,
      min_stream_order = min_stream_order,
      min_error = min_error,
      subset_label = subset_labels[q],
      output_filepath = trib_filepath
    )
    
    ################################################################################
    # STORE RESULTS IF REQUESTED
    ################################################################################
    
    if (return_values) {
      all_results[[q]] <- list(
        subset = current_subset,
        label = subset_labels[q],
        basin_assign_rescale = basin_assign_rescale,
        basin_assign_norm = basin_assign_norm,
        mgmt_result = mgmt_result  # Add management result
      )
    }
    
    message(paste("Completed processing quartile", q))
  }
  
  ################################################################################
  # RETURN RESULTS
  ################################################################################
  
  message(paste("=== DOY Quartile Analysis completed for", year, watershed, "==="))
  
  if (return_values) {
    return(all_results)
  } else {
    return(invisible(NULL))
  }
}

################################################################################
# SIMPLIFIED MANAGEMENT RIVER MAP FUNCTION
################################################################################

#' Create a simple management river map showing production proportion
create_simple_mgmt_map <- function(mgmt_result, edges, basin, gg_hist, year, watershed, 
                                   sensitivity_threshold, min_stream_order, subset_label, 
                                   output_filepath) {
  
  png(file = output_filepath, width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Set up plotting layout  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2, 
                                             heights = unit(c(0.7, 0.3), "npc"),
                                             widths = unit(c(0.6, 0.4), "npc"))))
  
  # Create color scale based on production proportion (like HUC maps)
  # Add production proportion to edges based on their management river
  edges$production_proportion <- NA
  for (i in 1:nrow(mgmt_result)) {
    mgmt_name <- mgmt_result$mgmt_river[i]
    prod_prop <- mgmt_result$production_proportion[i]
    edges$production_proportion[edges$mgmt_river == mgmt_name] <- prod_prop
  }
  
  # Filter to only managed edges for the map
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  # Convert to sf objects for ggplot
  managed_edges_sf <- st_as_sf(managed_edges)
  basin_sf <- st_as_sf(basin)
  
  # Main map plot showing management rivers colored by production proportion
  main_plot <- ggplot() +
    geom_sf(data = basin_sf, fill = "gray90", color = "gray70", size = 0.5) +
    geom_sf(data = managed_edges_sf, aes(color = production_proportion), size = 1.2) +
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
    coord_sf(datum = NA) +
    labs(
      title = paste0(subset_label, ": Management Rivers - ", watershed, " Watershed"),
      subtitle = paste("Year", year, "- Sensitivity:", sensitivity_threshold, 
                       "- Min Stream Order:", min_stream_order)
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
  
  # Create bar chart showing management rivers by production proportion
  # Apply watershed ordering (reverse for coord_flip so upstream appears at top)
  
  # First, ensure we have the ordering functions available
  if(!exists("get_watershed_order")) {
    get_watershed_order <- function() {
      # Updated to match the actual names in your data
      c("N. Fork Kusko", "E. Fork Kuskokwim", "S. Fork Kusko",
        "Upper Kusko Main", "Big River", "Takotna and Nixon Fork",
        "Tatlawiksuk", "Swift", "Stony", "Holitna and Hoholitna",
        "Middle Kusko Main", "George", "Oskakawlik", "Holokuk",
        "Aniak", "Tuluksak", "Kisaralik", "Kwethluk", "Johnson", "Lower Kusko")
    }
    
    apply_watershed_order <- function(data, mgmt_col = "mgmt_river", reverse_for_plots = FALSE) {
      standard_order <- get_watershed_order()
      if (reverse_for_plots) standard_order <- rev(standard_order)
      
      units_in_data <- unique(data[[mgmt_col]])
      final_order <- standard_order[standard_order %in% units_in_data]
      
      missing_units <- setdiff(units_in_data, standard_order)
      if (length(missing_units) > 0) {
        warning("Found management units not in standard order: ", paste(missing_units, collapse = ", "))
        cat("Units in data but not in standard order:\n")
        for(unit in missing_units) cat("  -", unit, "\n")
        final_order <- c(final_order, missing_units)
      }
      
      cat("Final ordering for", mgmt_col, ":\n")
      for(i in 1:length(final_order)) cat("  ", i, ".", final_order[i], "\n")
      
      data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
      return(data)
    }
  }
  
  # Apply watershed ordering (reverse for coord_flip so upstream appears at top)
  mgmt_result_ordered <- apply_watershed_order(mgmt_result, "mgmt_river", reverse_for_plots = TRUE)
  
  bar_plot <- ggplot(mgmt_result_ordered, 
                     aes(x = mgmt_river, y = production_proportion)) +
    geom_col(aes(fill = production_proportion), alpha = 0.9) +
    scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), guide = "none") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(title = "Production by Management River",
         x = "", y = "Production Proportion") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(5, 10, 5, 5, "mm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Plot components
  print(main_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(bar_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  # Add histogram if provided
  if (!is.null(gg_hist)) {
    gg_hist_clean <- gg_hist + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    print(gg_hist_clean, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  }
  
  dev.off()
  
  message(paste("Created Management map:", basename(output_filepath)))
}