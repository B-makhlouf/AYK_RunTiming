# cumulative_doy_total_analysis.R
# Functions for cumulative DOY quartile analysis normalized to total production
# UPDATED: Removed RawProduction folder and maps

library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(grid)
library(gridExtra)

# Source the required utility files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

#' Perform Cumulative DOY Quartile Analysis with Total Production Normalization
#'
#' @param year Character or numeric representing the year
#' @param watershed Character: "Kusko" or "Yukon"
#' @param sensitivity_threshold Numeric threshold for assignment filtering
#' @param min_error Minimum error value to use
#' @param min_stream_order Minimum stream order to include
#' @param HUC HUC level (e.g., 8, 10)
#' @param return_values Whether to return the calculated values
#' @return If return_values is TRUE, a list with results; otherwise NULL
Cumulative_DOY_Total_Analysis <- function(year, watershed, sensitivity_threshold, min_error, 
                                          min_stream_order = 3, HUC = 8, 
                                          return_values = FALSE) {
  
  # Generate identifier for output files
  identifier <- paste(year, watershed, sep = "_")
  
  # Load spatial data
  spatial_data <- load_spatial_data(watershed, HUC, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  Huc <- spatial_data$Huc
  
  # Load natal origins data
  natal_data <- load_natal_data(year, watershed)
  
  # Divide data into DOY quartiles
  quartile_data <- divide_doy_quartiles(natal_data)
  quartile_subsets <- quartile_data$subsets
  subset_labels <- quartile_data$labels
  
  # Extract isoscape prediction and error values
  pid_iso <- edges$iso_pred
  pid_isose <- edges$isose_pred
  
  # Calculate error values
  error <- calculate_error(pid_isose, min_error)
  
  # Set up watershed-specific priors
  priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
  
  # Create output directories (REMOVED RawProduction subdirectory)
  dir.create(here("Basin Maps/DOY_Cumulative_Total/HUC"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("Basin Maps/DOY_Cumulative_Total/Tribs"), showWarnings = FALSE, recursive = TRUE)
  
  # Create cumulative quartile subsets
  cumulative_subsets <- list()
  cumulative_labels <- list()
  
  # First quartile is just Q1
  cumulative_subsets[[1]] <- quartile_subsets[[1]]
  cumulative_labels[[1]] <- paste0("Cumulative DOY: ", subset_labels[1])
  
  # For remaining quartiles, add all previous data
  for (q in 2:length(quartile_subsets)) {
    # Combine data from quartiles 1 to q
    cumulative_subsets[[q]] <- do.call(rbind, quartile_subsets[1:q])
    
    # Create label that shows cumulative range
    q1_min <- min(quartile_subsets[[1]]$DOY, na.rm = TRUE)
    current_max <- max(quartile_subsets[[q]]$DOY, na.rm = TRUE)
    cumulative_labels[[q]] <- paste0("Cumulative DOY: ", ceiling(q1_min), "-", floor(current_max))
  }
  
  # First, calculate the TOTAL production for the full run (all data points)
  message("Calculating total production across all quartiles...")
  all_data <- do.call(rbind, quartile_subsets)
  all_assignment_matrix <- perform_assignment(
    all_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
  )
  
  # Calculate total basin production
  total_basin_assign_sum <- apply(all_assignment_matrix, 1, sum, na.rm = TRUE)
  grand_total_production <- sum(total_basin_assign_sum, na.rm = TRUE)
  
  # Return values storage
  if (return_values) {
    all_results <- list()
  }
  
  # Process each cumulative subset
  for (q in 1:length(cumulative_subsets)) {
    current_subset <- cumulative_subsets[[q]]
    
    # Skip empty subsets
    if (nrow(current_subset) == 0) {
      message(paste("Skipping", cumulative_labels[q], "because it contains no data"))
      next
    }
    
    # Create unique ID for this subset
    subset_id <- paste0(watershed, "_", year, "_CumulativeDOY_Q", q, "_Total")
    message(paste("Processing", cumulative_labels[[q]], "with", nrow(current_subset), "data points (as % of total run)"))
    
    # Perform assignment for this cumulative subset
    assignment_matrix <- perform_assignment(
      current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
    )
    
    # Calculate basin values for this cumulative subset
    cumulative_basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
    
    # Calculate DIRECT percentage of total run (not normalized within subset)
    cumulative_pct_of_total <- cumulative_basin_assign_sum / grand_total_production * 100
    
    # Process HUC data
    final_result <- process_huc_data(edges, basin, Huc, cumulative_basin_assign_sum, HUC)
    
    # Calculate percentages for each HUC
    final_result$total_production_all_quartiles <- grand_total_production
    final_result$percent_of_total_run <- (final_result$total_production / grand_total_production) * 100
    
    # Create improved histogram
    gg_hist <- create_doy_histogram(natal_data, current_subset, cumulative_labels[[q]])
    
    # Create HUC map with production per km
    huc_filepath <- file.path(here("Basin Maps/DOY_Cumulative_Total/HUC"), 
                              paste0(subset_id, "_HUC", HUC, ".png"))
    
    # Create HUC map
    png(file = huc_filepath, width = 12, height = 10, units = "in", res = 300, bg = "white")
    
    # Set up the plotting layout
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2, 
                                               heights = unit(c(0.7, 0.3), "npc"),
                                               widths = unit(c(0.6, 0.4), "npc"))))
    
    # Create main map plot with percent of total run
    main_plot <- ggplot() +
      geom_sf(data = final_result, aes(fill = percent_of_total_run), color = "white", size = 0.1) +
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
        title = paste0(cumulative_labels[[q]], ": Percent of Total Run - ", watershed, " Watershed"),
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
    
    # Sort HUCs by percent_of_total_run for the bar graph
    final_result_sorted <- final_result %>%
      arrange(desc(percent_of_total_run))
    
    # Create bar plot of percent of total run by HUC
    bargraph_percent <- ggplot(final_result_sorted, 
                               aes(x = reorder(Name, percent_of_total_run), y = percent_of_total_run)) +
      geom_col(aes(fill = percent_of_total_run), alpha = 0.9) +
      # Use YlOrRd color scale
      scale_fill_gradientn(
        colors = (brewer.pal(9, "YlOrRd")),
        name = "% of Total Run",
        labels = function(x) paste0(round(x, 1), "%")
      ) +
      coord_flip() +
      scale_y_continuous(
        expand = c(0, 0),
        labels = function(x) paste0(round(x, 1), "%")
      ) +
      labs(title = paste("Percent of Total Run by", "HUC", HUC),
           x = "",
           y = "Percent of Total Run") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(5, 10, 5, 5, "mm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    # If provided, make sure the gg_hist has white background
    if (!is.null(gg_hist)) {
      gg_hist <- enforce_histogram_limits(gg_hist) + 
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        )
    }
    
    # Plot main map in upper left
    print(main_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    
    # Plot production bar chart in upper right
    print(bargraph_percent, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    
    # If we have a histogram to include, add it at the bottom spanning both columns
    if (!is.null(gg_hist)) {
      print(gg_hist, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    }
    
    dev.off()
    
    # REMOVED: Create HUC map with raw production values
    # The entire section that created raw_huc_filepath and raw production maps has been removed
    
    # Create tributary map
    trib_filepath <- file.path(here("Basin Maps/DOY_Cumulative_Total/Tribs"), 
                               paste0(subset_id, ".png"))
    
    # Open PNG for tributary map
    png(file = trib_filepath, width = 9, height = 8, units = "in", res = 300, bg = "white")
    
    # Use the YlOrRd palette with 9 colors expanded to 10
    pallete <- brewer.pal(9, "YlOrRd")
    pallete_expanded <- colorRampPalette(pallete)(10)
    
    # Instead of normalizing values, use actual percentages of total for coloring
    # Find the maximum percentage for scaling
    max_pct <- max(cumulative_pct_of_total, na.rm = TRUE)
    
    # Scale percentages for coloring - use breakpoints based on percentile of actual values
    colcode <- rep("gray60", length(cumulative_pct_of_total))
    colcode[cumulative_pct_of_total == 0] <- 'white'
    
    # Use data-driven breakpoints for better visualization
    # Example: 0-10% of max, 10-20% of max, etc.
    for (i in 1:10) {
      lower <- (i-1) * max_pct / 10
      upper <- i * max_pct / 10
      colcode[cumulative_pct_of_total > lower & cumulative_pct_of_total <= upper] <- pallete_expanded[i]
    }
    
    # Override with grays for areas we don't want to show
    colcode[which(priors$StreamOrderPrior == 0)] <- 'gray60'
    colcode[which(priors$pid_prior == 0)] <- 'gray60'
    
    # Set linewidths based on stream order
    if (watershed == "Yukon") {
      stream_order_lwd <- edges$Str_Order
      linewidths <- rep(1, length(stream_order_lwd))
      linewidths <- ifelse(stream_order_lwd == 9, 3.7, linewidths)
      linewidths <- ifelse(stream_order_lwd == 8, 2.5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 7, 1.7, linewidths)
      linewidths <- ifelse(stream_order_lwd == 6, 1.5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 5, 1, linewidths)
      linewidths <- ifelse(stream_order_lwd == 4, 1, linewidths)
      linewidths <- ifelse(stream_order_lwd == 3, 1, linewidths)
    } else {
      stream_order_lwd <- edges$Str_Order
      linewidths <- rep(1, length(stream_order_lwd))
      linewidths <- ifelse(stream_order_lwd == 9, 5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 8, 4, linewidths)
      linewidths <- ifelse(stream_order_lwd == 7, 3, linewidths)
      linewidths <- ifelse(stream_order_lwd == 6, 2, linewidths)
      linewidths <- ifelse(stream_order_lwd == 5, 1.8, linewidths)
      linewidths <- ifelse(stream_order_lwd == 4, 1.5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 3, 1, linewidths)
    }
    
    # Generate title
    plot_title <- paste0(
      cumulative_labels[[q]], " (% of Total Run)\n",
      "Year:", year, 
      " River:", watershed,
      "\nThreshold:", sensitivity_threshold, 
      " Min Stream Order:", min_stream_order,
      " Min Error:", min_error
    )
    
    # Adjust plot margins
    par(mar = c(8, 4, 4, 2), bg = "white")
    
    # Plot the basin and edges
    plot(st_geometry(basin), col = 'gray60', border = 'gray60', main = plot_title, bg = "white")
    plot(st_geometry(edges), col = colcode, pch = 16, axes = FALSE, add = TRUE, lwd = linewidths)
    
    # Create legend breaks based on actual percentage values
    legend_breaks <- seq(0, max_pct, length.out = 11)
    legend_labels <- sapply(1:10, function(i) {
      sprintf("%.1f-%.1f%%", legend_breaks[i], legend_breaks[i+1])
    })
    
    # Add legend
    legend("topleft", 
           legend = legend_labels, 
           col = pallete_expanded, 
           lwd = 5, 
           title = "Percent of Total Run", 
           bty = "n",
           bg = "white")
    
    # Add histogram to the plot if provided
    if (!is.null(gg_hist)) {
      # Modify the histogram specifically for grid viewport use
      limited_hist <- gg_hist +
        scale_x_continuous(limits = c(140, 200)) +
        scale_y_continuous(limits = c(0, 0.1)) +
        coord_cartesian(xlim = c(140, 200), ylim = c(0, 0.1), expand = FALSE) +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(0, 0, 0, 0)
        )
      
      # Create viewport with explicit scaling
      vp_hist <- viewport(
        x = 0.5, y = 0.05, 
        width = 0.7, height = 0.2, 
        just = c("center", "bottom")
      )
      
      # Print the modified histogram
      print(limited_hist, vp = vp_hist)
    }
    
    dev.off()
    
    # Reset par to default
    par(mar = c(5, 4, 4, 2) + 0.1, bg = "white")
    
    # Store results if needed
    if (return_values) {
      all_results[[q]] <- list(
        subset = current_subset,
        label = cumulative_labels[[q]],
        percent_of_total_run = cumulative_pct_of_total,
        total_basin_assign_sum = total_basin_assign_sum,
        grand_total_production = grand_total_production,
        huc_result = final_result
      )
    }
  }
  
  if (return_values) {
    return(all_results)
  } else {
    return(invisible(NULL))
  }
}