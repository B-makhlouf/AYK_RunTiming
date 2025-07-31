# visualization.R
# Visualization utilities for the DOY timing analysis with PNG output
# NOTE: HUC map functions removed - only management river and tributary functions remain

library(ggplot2)
library(RColorBrewer)
library(scales)
library(grid)
library(gridExtra)
library(viridis)

#' Create a DOY distribution histogram highlighting a specific subset
#'
#' @param full_dataset Full dataset containing DOY and dailyCPUEprop
#' @param current_subset Subset of data to highlight
#' @param title Optional custom title for the plot
#' @return ggplot object with the histogram
create_doy_histogram <- function(full_dataset, current_subset, title = NULL) {
  # Set default title if not provided
  if (is.null(title)) {
    title <- "DOY Distribution"
  }
  
  # Create a custom color for highlighting the subset
  highlight_color <- "tomato"
  background_color <- "gray70"
  
  # Create a function to convert DOY to date for the x-axis labels
  doy_to_date <- function(doy, year = 2024) {  # Use 2024 as default (leap year)
    as.Date(doy - 1, origin = paste0(year, "-01-01"))
  }
  
  # Create custom x-axis breaks and labels
  doy_breaks <- seq(140, 210, by = 10)
  
  # Create the histogram with FIXED axis limits and coordinate system
  ggplot() + 
    # First plot the full dataset as background with reduced opacity
    # Always show the complete CPUE curve regardless of subset
    geom_line(data = full_dataset, aes(x = DOY, y = dailyCPUEprop), 
              color = "gray40", linewidth = 1, alpha = 0.5) +
    geom_ribbon(data = full_dataset, aes(x = DOY, ymin = 0, ymax = dailyCPUEprop), 
                fill = background_color, alpha = 0.3) +
    
    # Then overlay the current subset with higher opacity and different color
    geom_line(data = current_subset, aes(x = DOY, y = dailyCPUEprop), 
              color = "black", linewidth = 2) +
    geom_ribbon(data = current_subset, aes(x = DOY, ymin = 0, ymax = dailyCPUEprop), 
                fill = highlight_color, alpha = 0.7) +
    
    # Add vertical lines to show the subset boundaries
    geom_vline(xintercept = c(
      min(current_subset$DOY), 
      max(current_subset$DOY)
    ), linetype = "dashed", color = "darkred") +
    
    # IMPORTANT: Set consistent axis ranges with custom breaks and labels
    scale_x_continuous(limits = c(140, 210),
                       breaks = doy_breaks,
                       labels = function(x) {
                         # Create two-line labels with DOY and date
                         paste0(x, "\n", format(doy_to_date(x), "%b %d"))
                       }) +
    scale_y_continuous(limits = c(0, 0.1)) +
    # Force the coordinate system to respect our limits without expansion
    coord_cartesian(xlim = c(140, 210), ylim = c(0, 0.1), expand = FALSE) +
    
    # Add labels and theme
    labs(
      title = title,
      subtitle = "Current subset highlighted in red",
      x = "Day of Year (Date)", 
      y = "Daily CPUE Proportion"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

#' Enforce consistent histogram limits
#'
#' @param gg_hist ggplot object to modify
#' @return Modified ggplot object with enforced limits
enforce_histogram_limits <- function(gg_hist) {
  if (!is.null(gg_hist)) {
    # Create a function to convert DOY to date for the x-axis labels
    doy_to_date <- function(doy, year = 2024) {  # Use 2024 as default (leap year)
      as.Date(doy - 1, origin = paste0(year, "-01-01"))
    }
    
    # Create custom x-axis breaks and labels
    doy_breaks <- seq(140, 210, by = 10)
    
    gg_hist <- gg_hist + 
      scale_x_continuous(limits = c(140, 210), 
                         breaks = doy_breaks,
                         labels = function(x) {
                           # Create two-line labels with DOY and date
                           paste0(x, "\n", format(doy_to_date(x), "%b %d"))
                         },
                         expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 0.1), 
                         breaks = seq(0, 0.1, by = 0.02),
                         expand = c(0, 0)) +
      coord_cartesian(xlim = c(140, 210), ylim = c(0, 0.1), expand = FALSE) +
      labs(x = "Day of Year (Date)") +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(hjust = 1)
      )
  }
  return(gg_hist)
}

#' Create and save a DOY tributary map as PNG
#'
#' @param basin SF object with basin boundary
#' @param edges SF object with stream edges
#' @param basin_assign_norm Vector of normalized basin assignments
#' @param StreamOrderPrior Vector of stream order priors
#' @param pid_prior Vector of watershed priors
#' @param gg_hist ggplot object with histogram to overlay
#' @param year Year for the title
#' @param watershed Watershed name for the title
#' @param sensitivity_threshold Sensitivity threshold used
#' @param min_stream_order Minimum stream order used
#' @param min_error Minimum error used
#' @param subset_label Optional label for the subset
#' @param output_filepath Path to save the PNG
#' @return Invisibly returns the output filepath
create_doy_tributary_map <- function(basin, edges, basin_assign_norm, StreamOrderPrior, 
                                     pid_prior, gg_hist, year, watershed, 
                                     sensitivity_threshold, min_stream_order, min_error,
                                     subset_label = NULL, output_filepath) {
  
  # Update file extension from pdf to png if needed
  output_filepath <- sub("\\.pdf$", ".png", output_filepath)
  
  # Open PNG for tributary map - higher resolution
  png(file = output_filepath, width = 9, height = 8, units = "in", res = 300, bg = "white")
  
  # Use the YlOrRd palette with 9 colors expanded to 10
  pallete <- brewer.pal(9, "YlOrRd")
  pallete_expanded <- colorRampPalette(pallete)(10)
  
  # Color coding with bins at every 0.1
  colcode <- rep("gray60", length(basin_assign_norm))
  colcode[basin_assign_norm == 0] <- 'white'
  colcode[basin_assign_norm > 0 & basin_assign_norm <= 0.2] <- pallete_expanded[1]
  colcode[basin_assign_norm > 0.2 & basin_assign_norm <= 0.4] <- pallete_expanded[4]
  colcode[basin_assign_norm > 0.4 & basin_assign_norm <= 0.6] <- pallete_expanded[5]
  colcode[basin_assign_norm > 0.6 & basin_assign_norm <= 0.7] <- pallete_expanded[7]
  colcode[basin_assign_norm > 0.7 & basin_assign_norm <= 0.8] <- pallete_expanded[8]
  colcode[basin_assign_norm > 0.8 & basin_assign_norm <= 0.9] <- pallete_expanded[9]
  colcode[basin_assign_norm > 0.9 & basin_assign_norm <= 1.0] <- pallete_expanded[10]
  colcode[which(StreamOrderPrior == 0)] <- 'gray60'
  colcode[which(pid_prior == 0)] <- 'gray60'
  
  if (watershed == "Yukon") {
    # Set linewidths based on stream order and probability
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
    # Set linewidths based on stream order and probability
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
    ifelse(is.null(subset_label), "", paste0(subset_label, "\n")),
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
  
  # Add legend
  legend("topleft", 
         legend = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                    "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"), 
         col = pallete_expanded, 
         lwd = 5, 
         title = "Relative posterior density", 
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
  
  message(paste("Created PNG tributary map:", output_filepath))
  invisible(output_filepath)
}