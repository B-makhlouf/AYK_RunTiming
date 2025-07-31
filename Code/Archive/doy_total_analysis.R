################################################################################
# DOY_TOTAL_ANALYSIS.R - DOY ANALYSIS WITH TOTAL RUN NORMALIZATION
################################################################################
# PURPOSE: Analyzes salmon run timing quartiles as proportion of TOTAL ANNUAL RUN
# OUTPUT: Shows what % of the ENTIRE year's run each quartile represents
# KEY DIFFERENCE: Normalizes to total run, not within quartile
################################################################################

library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)

# Source required utilities  
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

################################################################################
# MAIN DOY TOTAL ANALYSIS FUNCTION
################################################################################

#' Perform DOY Total Analysis - normalized to entire annual run
DOY_Total_Analysis <- function(year, watershed, sensitivity_threshold, min_error, 
                               min_stream_order = 3, 
                               return_values = FALSE) {
  
  message(paste("=== Starting DOY Total Analysis for", year, watershed, "==="))
  
  ################################################################################
  # SETUP: LOAD DATA AND PARAMETERS
  ################################################################################
  
  # Load spatial and natal data
  spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  
  natal_data <- load_natal_data(year, watershed)
  
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
  # CALCULATE TOTAL ANNUAL PRODUCTION (FOR NORMALIZATION)
  ################################################################################
  
  message("Calculating total production across all quartiles for normalization...")
  
  # Combine all quartiles to get total annual data  
  all_data <- do.call(rbind, quartile_subsets)
  
  # Perform assignment on entire annual dataset
  all_assignment_matrix <- perform_assignment(
    all_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
  )
  
  # Calculate total basin production for the entire year
  total_basin_assign_sum <- apply(all_assignment_matrix, 1, sum, na.rm = TRUE)
  grand_total_production <- sum(total_basin_assign_sum, na.rm = TRUE)
  
  message(paste("Grand total production:", round(grand_total_production, 6)))
  
  ################################################################################
  # CREATE OUTPUT DIRECTORIES
  ################################################################################
  
  dir.create(here("Basin Maps/DOY_Total/Tribs"), showWarnings = FALSE, recursive = TRUE)
  
  # Storage for return values
  if (return_values) {
    all_results <- list()
  }
  
  ################################################################################
  # PROCESS EACH QUARTILE (NORMALIZED TO TOTAL RUN)
  ################################################################################
  
  for (q in 1:length(quartile_subsets)) {
    current_subset <- quartile_subsets[[q]]
    
    # Skip empty quartiles
    if (nrow(current_subset) == 0) {
      message(paste("Skipping Q", q, "because it contains no data"))
      next
    }
    
    message(paste("Processing DOY Q", q, "with", nrow(current_subset), "data points (as % of total run)"))
    
    ################################################################################
    # PERFORM ASSIGNMENT FOR THIS QUARTILE
    ################################################################################
    
    subset_id <- paste0(watershed, "_", year, "_DOY_Q", q, "_Total")
    
    # Perform assignment for this quartile
    assignment_matrix <- perform_assignment(
      current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Calculate basin values for this quartile
    quartile_basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
    
    ################################################################################
    # CALCULATE PERCENTAGE OF TOTAL RUN (KEY DIFFERENCE FROM REGULAR DOY)
    ################################################################################
    
    # Calculate as DIRECT proportion of total run (not normalized within quartile)
    quartile_pct_of_total <- quartile_basin_assign_sum / grand_total_production * 100
    
    ################################################################################
    # CREATE VISUALIZATION PLOTS
    ################################################################################
    
    # Create DOY histogram
    gg_hist <- create_doy_histogram(natal_data, current_subset, subset_labels[q])
    
    ################################################################################
    # CREATE TRIBUTARY MAP WITH TOTAL RUN PERCENTAGES  
    ################################################################################
    
    trib_filepath <- file.path(here("Basin Maps/DOY_Total/Tribs"), 
                               paste0(subset_id, ".png"))
    
    create_total_tributary_map(basin, edges, quartile_pct_of_total, priors, 
                               gg_hist, year, watershed, sensitivity_threshold, 
                               min_stream_order, min_error, subset_labels[q], 
                               trib_filepath)
    
    ################################################################################
    # STORE RESULTS IF REQUESTED
    ################################################################################
    
    if (return_values) {
      all_results[[q]] <- list(
        subset = current_subset,
        label = subset_labels[q],
        percent_of_total_run = quartile_pct_of_total,
        total_basin_assign_sum = total_basin_assign_sum,
        grand_total_production = grand_total_production
      )
    }
  }
  
  ################################################################################
  # RETURN RESULTS
  ################################################################################
  
  message(paste("=== DOY Total Analysis completed for", year, watershed, "==="))
  
  if (return_values) {
    return(all_results)
  } else {
    return(invisible(NULL))
  }
}

################################################################################
# HELPER FUNCTIONS FOR TOTAL RUN MAPS
################################################################################

#' Create tributary map showing percent of total run  
create_total_tributary_map <- function(basin, edges, quartile_pct_of_total, priors, 
                                       gg_hist, year, watershed, sensitivity_threshold, 
                                       min_stream_order, min_error, subset_label, 
                                       output_filepath) {
  
  png(file = output_filepath, width = 9, height = 8, units = "in", res = 300, bg = "white")
  
  # Color setup
  pallete <- brewer.pal(9, "YlOrRd")
  pallete_expanded <- colorRampPalette(pallete)(10)
  
  max_pct <- max(quartile_pct_of_total, na.rm = TRUE)
  
  # Color coding based on actual percentages
  colcode <- rep("gray60", length(quartile_pct_of_total))
  colcode[quartile_pct_of_total == 0] <- 'white'
  
  for (i in 1:10) {
    lower <- (i-1) * max_pct / 10
    upper <- i * max_pct / 10
    colcode[quartile_pct_of_total > lower & quartile_pct_of_total <= upper] <- pallete_expanded[i]
  }
  
  # Override with grays for excluded areas
  colcode[which(priors$StreamOrderPrior == 0)] <- 'gray60'
  colcode[which(priors$pid_prior == 0)] <- 'gray60'
  
  # Set line widths by stream order
  stream_order_lwd <- edges$Str_Order
  linewidths <- rep(1, length(stream_order_lwd))
  if (watershed == "Yukon") {
    linewidths <- ifelse(stream_order_lwd == 9, 3.7, linewidths)
    linewidths <- ifelse(stream_order_lwd == 8, 2.5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 7, 1.7, linewidths)
    linewidths <- ifelse(stream_order_lwd == 6, 1.5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 5, 1, linewidths)
  } else {
    linewidths <- ifelse(stream_order_lwd == 9, 5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 8, 4, linewidths)
    linewidths <- ifelse(stream_order_lwd == 7, 3, linewidths)
    linewidths <- ifelse(stream_order_lwd == 6, 2, linewidths)
    linewidths <- ifelse(stream_order_lwd == 5, 1.8, linewidths)
    linewidths <- ifelse(stream_order_lwd == 4, 1.5, linewidths)
  }
  
  # Create plot
  plot_title <- paste0(subset_label, " (% of Total Run)\n",
                       "Year:", year, " River:", watershed,
                       "\nThreshold:", sensitivity_threshold, 
                       " Min Stream Order:", min_stream_order, " Min Error:", min_error)
  
  par(mar = c(8, 4, 4, 2), bg = "white")
  plot(st_geometry(basin), col = 'gray60', border = 'gray60', main = plot_title, bg = "white")
  plot(st_geometry(edges), col = colcode, pch = 16, axes = FALSE, add = TRUE, lwd = linewidths)
  
  # Legend
  legend_breaks <- seq(0, max_pct, length.out = 11)
  legend_labels <- sapply(1:10, function(i) {
    sprintf("%.1f-%.1f%%", legend_breaks[i], legend_breaks[i+1])
  })
  
  legend("topleft", legend = legend_labels, col = pallete_expanded, lwd = 5, 
         title = "Percent of Total Run", bty = "n", bg = "white")
  
  # Add histogram if provided
  if (!is.null(gg_hist)) {
    limited_hist <- gg_hist + scale_x_continuous(limits = c(140, 200)) +
      scale_y_continuous(limits = c(0, 0.1)) + coord_cartesian(xlim = c(140, 200), ylim = c(0, 0.1), expand = FALSE) +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA), plot.margin = margin(0, 0, 0, 0))
    
    vp_hist <- viewport(x = 0.5, y = 0.05, width = 0.7, height = 0.2, just = c("center", "bottom"))
    print(limited_hist, vp = vp_hist)
  }
  
  dev.off()
  par(mar = c(5, 4, 4, 2) + 0.1, bg = "white")
}