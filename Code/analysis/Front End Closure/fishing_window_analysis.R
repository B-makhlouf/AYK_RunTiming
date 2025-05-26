# fishing_window_analysis.R
# Analysis of salmon production within specific fishing window (June 1-11)

library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(grid)
library(gridExtra)

# Source the required utility files
source(here("code/utils/spatial_utils.R"))
source(here("code/utils/visualization.R"))
source(here("code/assignment.R"))

#' Analyze salmon production within a fishing window
#'
#' @param year Character or numeric representing the year
#' @param watershed Character: "Kusko" or "Yukon"
#' @param sensitivity_threshold Numeric threshold for assignment filtering
#' @param min_error Minimum error value to use
#' @param min_stream_order Minimum stream order to include
#' @param HUC HUC level (e.g., 8, 10)
#' @param start_date Start date for the fishing window
#' @param end_date End date for the fishing window
#' @param return_values Whether to return the calculated values
#' @return If return_values is TRUE, a list with results; otherwise NULL
analyze_fishing_window <- function(year, watershed, sensitivity_threshold, min_error, 
                                   min_stream_order = 3, HUC = 8,
                                   start_date = "06-01", end_date = "06-11",
                                   return_values = FALSE) {
  
  # Convert dates to DOY (Day of Year)
  # Using 2020 as a reference year (which was a leap year) for consistent DOY calculations
  start_doy <- as.numeric(format(as.Date(paste0("2020-", start_date)), "%j"))
  end_doy <- as.numeric(format(as.Date(paste0("2020-", end_date)), "%j"))
  
  # Print the DOY values for verification
  message(paste("Analyzing salmon run within fishing window (DOY", start_doy, "-", end_doy, ")"))
  
  # Generate identifier for output files
  identifier <- paste(year, watershed, sep = "_")
  
  # Load spatial data
  spatial_data <- load_spatial_data(watershed, HUC, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  Huc <- spatial_data$Huc
  
  # Load natal origins data
  natal_data <- load_natal_data(year, watershed)
  
  # Extract isoscape prediction and error values
  pid_iso <- edges$iso_pred
  pid_isose <- edges$isose_pred
  
  # Calculate error values
  error <- calculate_error(pid_isose, min_error)
  
  # Set up watershed-specific priors
  priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
  
  # Create output directories
  dir.create(here("Basin Maps/Fishing_Window/HUC"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("Basin Maps/Fishing_Window/HUC/RawProduction"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("Basin Maps/Fishing_Window/Tribs"), showWarnings = FALSE, recursive = TRUE)
  
  # Create window-specific subset of data (fish that arrive within the fishing window)
  window_data <- natal_data %>% filter(DOY >= start_doy & DOY <= end_doy)
  
  # Skip if no data in the window period
  if (nrow(window_data) == 0) {
    message(paste("No data found within the fishing window (DOY", start_doy, "-", end_doy, ") for", year, watershed))
    return(NULL)
  }
  
  # Calculate total run statistics (all data points)
  message("Calculating total production across entire run...")
  all_assignment_matrix <- perform_assignment(
    natal_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
  )
  
  # Calculate total basin production
  total_basin_assign_sum <- apply(all_assignment_matrix, 1, sum, na.rm = TRUE)
  grand_total_production <- sum(total_basin_assign_sum, na.rm = TRUE)
  
  # Calculate window statistics
  message(paste("Calculating production within fishing window (DOY", start_doy, "-", end_doy, ")..."))
  window_assignment_matrix <- perform_assignment(
    window_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold
  )
  
  # Calculate window basin values
  window_basin_assign_sum <- apply(window_assignment_matrix, 1, sum, na.rm = TRUE)
  window_total_production <- sum(window_basin_assign_sum, na.rm = TRUE)
  
  # Calculate DIRECT proportion of total run production (not normalized within window)
  # This ensures each value directly represents what percent of the total run it contributes
  window_pct_of_total <- window_basin_assign_sum / grand_total_production * 100
  
  # Calculate overall percentage of run that falls within the fishing window
  overall_window_pct <- window_total_production / grand_total_production * 100
  
  # Print summary statistics
  message(paste("Overall percentage of run within fishing window (DOY", start_doy, "-", end_doy, "):", 
                round(overall_window_pct, 2), "%"))
  
  # Create label for the analysis
  window_label <- paste0("Fishing Window (DOY ", start_doy, "-", end_doy, ", ", 
                         format(as.Date(paste0("2020-", start_date)), "%b %d"), "-",
                         format(as.Date(paste0("2020-", end_date)), "%b %d"), ")")
  subset_id <- paste0(watershed, "_", year, "_FishingWindow_", start_doy, "_", end_doy)
  
  # Process HUC data for full run to calculate total percentages
  all_hucs_result <- process_huc_data(edges, basin, Huc, total_basin_assign_sum, HUC)
  all_hucs_result$percent_of_total_run <- (all_hucs_result$total_production / grand_total_production) * 100
  
  # Process HUC data for window period
  window_result <- process_huc_data(edges, basin, Huc, window_basin_assign_sum, HUC)
  window_result$total_production_all_data <- grand_total_production
  window_result$percent_of_total_run <- (window_result$total_production / grand_total_production) * 100
  
  # Create improved histogram highlighting the fishing window
  gg_hist <- create_doy_histogram(natal_data, window_data, window_label)
  
  # Create HUC map with production per km
  huc_filepath <- file.path(here("Basin Maps/Fishing_Window/HUC"), 
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
    geom_sf(data = window_result, aes(fill = percent_of_total_run), color = "white", size = 0.1) +
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
                       round(overall_window_pct, 1), "%")
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
  
  # Create the updated histogram showing total vs. window contribution
  huc_histogram <- create_huc_total_quartile_histogram(
    final_result = window_result,
    all_hucs_data = all_hucs_result,
    quartile_label = window_label
  )
  
  # Print the plots to the PNG
  print(main_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(huc_histogram, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  # Add the DOY histogram at the bottom spanning both columns
  if (!is.null(gg_hist)) {
    gg_hist <- enforce_histogram_limits(gg_hist) + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      # Add vertical lines at the window boundaries
      geom_vline(xintercept = start_doy, color = "blue", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = end_doy, color = "blue", linetype = "dashed", linewidth = 1) +
      annotate("text", x = start_doy - 2, y = 0.09, 
               label = paste0("Start: DOY ", start_doy, "\n", format(as.Date(paste0("2020-", start_date)), "%b %d")),
               color = "blue", hjust = 1) +
      annotate("text", x = end_doy + 2, y = 0.09, 
               label = paste0("End: DOY ", end_doy, "\n", format(as.Date(paste0("2020-", end_date)), "%b %d")),
               color = "blue", hjust = 0)
    
    print(gg_hist, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  }
  
  dev.off()
  
  # Create tributary map
  trib_filepath <- file.path(here("Basin Maps/Fishing_Window/Tribs"), 
                             paste0(subset_id, ".png"))
  
  # Open PNG for tributary map
  png(file = trib_filepath, width = 9, height = 8, units = "in", res = 300, bg = "white")
  
  # Use the YlOrRd palette with 9 colors expanded to 10
  pallete <- brewer.pal(9, "YlOrRd")
  pallete_expanded <- colorRampPalette(pallete)(10)
  
  # Instead of normalizing values, use actual percentages of total for coloring
  # Find the maximum percentage for scaling
  max_pct <- max(window_pct_of_total, na.rm = TRUE)
  
  # Scale percentages for coloring - use breakpoints based on percentile of actual values
  colcode <- rep("gray60", length(window_pct_of_total))
  colcode[window_pct_of_total == 0] <- 'white'
  
  # Use data-driven breakpoints for better visualization
  # Example: 0-10% of max, 10-20% of max, etc.
  for (i in 1:10) {
    lower <- (i-1) * max_pct / 10
    upper <- i * max_pct / 10
    colcode[window_pct_of_total > lower & window_pct_of_total <= upper] <- pallete_expanded[i]
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
    window_label, " (% of Total Run)\n",
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
      # Add vertical lines at the fishing window boundaries
      geom_vline(xintercept = start_doy, color = "blue", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = end_doy, color = "blue", linetype = "dashed", linewidth = 1) +
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
  
  # Create a summary stats data frame
  summary_stats <- data.frame(
    year = year,
    watershed = watershed,
    window_start_doy = start_doy,
    window_end_doy = end_doy,
    window_start_date = format(as.Date(paste0("2020-", start_date)), "%b %d"),
    window_end_date = format(as.Date(paste0("2020-", end_date)), "%b %d"),
    total_production = grand_total_production,
    window_production = window_total_production,
    window_pct = overall_window_pct
  )
  
  # Return values if requested
  if (return_values) {
    return(list(
      window_data = window_data,
      label = window_label,
      percent_of_total_run = window_pct_of_total,
      total_basin_assign_sum = total_basin_assign_sum,
      grand_total_production = grand_total_production,
      window_basin_assign_sum = window_basin_assign_sum,
      window_total_production = window_total_production,
      overall_window_pct = overall_window_pct,
      huc_result = window_result,
      all_hucs_result = all_hucs_result,
      summary_stats = summary_stats
    ))
  } else {
    return(invisible(NULL))
  }
}

#' Run the fishing window analysis for multiple years (Kuskokwim watershed only)
#'
#' @param years Vector of years to process
#' @param start_date Start date for the fishing window
#' @param end_date End date for the fishing window
#' @param export_results Whether to export results to CSV
#' @return Data frame with summary statistics
run_fishing_window_analysis <- function(years, 
                                        start_date = "06-01", end_date = "06-11",
                                        export_results = TRUE) {
  # Create a list to store results
  all_results <- list()
  summary_stats_list <- list()
  
  # Create output directories
  results_dir <- here("Analysis_Results/Fishing_Window")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Set Kuskokwim-specific parameters
  watershed <- "Kusko"
  sensitivity_threshold <- 0.7
  min_error <- 0.0006
  min_stream_order <- 3
  
  # Process each year
  for (year in years) {
    message(paste("\nAnalyzing fishing window for Kuskokwim,", year))
    
    result <- analyze_fishing_window(
      year = year,
      watershed = watershed,
      sensitivity_threshold = sensitivity_threshold,
      min_error = min_error,
      min_stream_order = min_stream_order,
      HUC = 8,
      start_date = start_date,
      end_date = end_date,
      return_values = TRUE
    )
    
    if (!is.null(result)) {
      all_results[[paste(year, watershed, sep = "_")]] <- result
      summary_stats_list[[paste(year, watershed, sep = "_")]] <- result$summary_stats
    }
  }
  
  # Combine all summary statistics
  if (length(summary_stats_list) > 0) {
    all_summary_stats <- bind_rows(summary_stats_list)
    
    # Export results if requested
    if (export_results) {
      # Export summary statistics
      summary_filepath <- file.path(results_dir, "fishing_window_summary.csv")
      write.csv(all_summary_stats, summary_filepath, row.names = FALSE)
      message(paste("Saved summary statistics to:", summary_filepath))
      
      # Export HUC data if available
      if (length(all_results) > 0) {
        # Combine HUC data from all results
        all_huc_data <- bind_rows(lapply(names(all_results), function(name) {
          result <- all_results[[name]]
          huc_data <- st_drop_geometry(result$huc_result)
          parts <- strsplit(name, "_")[[1]]
          huc_data$year <- parts[1]
          return(huc_data)
        }))
        
        # Save HUC data
        huc_filepath <- file.path(results_dir, "fishing_window_huc_data.csv")
        write.csv(all_huc_data, huc_filepath, row.names = FALSE)
        message(paste("Saved HUC data to:", huc_filepath))
      }
    }
    
    return(all_summary_stats)
  } else {
    warning("No results were generated")
    return(NULL)
  }
}

#' Create a summary plot of fishing window percentages across years for Kuskokwim
#'
#' @param summary_stats Data frame with summary statistics
#' @param output_dir Directory to save output file
#' @return Path to the saved PNG file
create_summary_plot <- function(summary_stats, output_dir = here("Basin Maps/Fishing_Window")) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Format the window dates for the title
  window_dates <- paste(unique(summary_stats$window_start_date), "-", unique(summary_stats$window_end_date))
  
  # Create a plot showing window percentages by year
  p <- ggplot(summary_stats, aes(x = as.factor(year), y = window_pct)) +
    geom_bar(stat = "identity", fill = "#E41A1C") +
    geom_text(aes(label = sprintf("%.1f%%", window_pct)), 
              vjust = -0.5, size = 4) +
    labs(
      title = paste("Percentage of Total Salmon Run During Fishing Window (", window_dates, ")"),
      subtitle = "Kuskokwim Watershed",
      x = "Year",
      y = "Percentage of Total Run (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Save the plot
  plot_filepath <- file.path(output_dir, "fishing_window_percentage_summary.png")
  ggsave(plot_filepath, p, width = 10, height = 6, dpi = 300)
  
  message(paste("Created summary plot:", plot_filepath))
  return(plot_filepath)
}

#' Create a HUC comparison plot across years for Kuskokwim
#'
#' @param all_huc_data Data frame with HUC data from all years
#' @param output_dir Directory to save output file
#' @return Path to the saved PNG file
create_huc_comparison_plot <- function(all_huc_data, output_dir = here("Basin Maps/Fishing_Window")) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Focus on the top contributing HUCs
  top_hucs <- all_huc_data %>%
    group_by(Name) %>%
    summarize(avg_pct = mean(percent_of_total_run, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(avg_pct)) %>%
    slice_head(n = 10) %>%
    pull(Name)
  
  # Filter data for the top HUCs
  top_huc_data <- all_huc_data %>%
    filter(Name %in% top_hucs) %>%
    # Rename year for better plotting
    rename(Year = year) %>% 
    # Order HUCs by average contribution
    mutate(Name = factor(Name, levels = top_hucs))
  
  # Create a plot showing HUC contributions by year
  p <- ggplot(top_huc_data, aes(x = as.factor(Year), y = percent_of_total_run, fill = Name)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set3", name = "HUC") +
    labs(
      title = "Top 10 HUC Contributions During Fishing Window",
      subtitle = "Percentage of Total Run by HUC and Year - Kuskokwim Watershed",
      x = "Year",
      y = "Percentage of Total Run (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10)
    )
  
  # Save the plot
  plot_filepath <- file.path(output_dir, "fishing_window_huc_comparison.png")
  ggsave(plot_filepath, p, width = 12, height = 8, dpi = 300)
  
  message(paste("Created HUC comparison plot:", plot_filepath))
  return(plot_filepath)
}

# Example usage - modified to ONLY generate maps, no summary figures
# Run the analysis for years 2017-2021
run_fishing_window_analysis(
  years = c("2017", "2018", "2019", "2020", "2021"),
  start_date = "06-01", 
  end_date = "06-11",
  export_results = TRUE  # Still export the data files, but don't create summary plots
)


# create_summary_plot(summary_stats)
# create_huc_comparison_plot(all_huc_data)