################################################################################
# COMPLETE CUMULATIVE DISTRIBUTION ANALYSIS SCRIPT - FIXED NAMES
################################################################################
# PURPOSE: Calculate cumulative distribution plots for each management unit
#          with data points every 3 days - READY TO RUN
# FIXED: Corrected management unit names to match actual data
################################################################################

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(RColorBrewer)
  library(tidyr)
  library(scales)
  library(gridExtra)
  library(grid)
})

# Source required utilities
cat("Loading spatial utilities and assignment functions...\n")
tryCatch({
  source(here("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Code/utils/spatial_utils.R"))
  source(here("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Code/assignment.R"))
}, error = function(e) {
  stop("Could not load required utility files. Please ensure spatial_utils.R and assignment.R are available.\nError: ", e$message)
})

################################################################################
# FUNCTION DEFINITIONS
################################################################################

#' Create CPUE histogram for a specific year with fixed June 1 - July 30 axis
create_cpue_histogram <- function(natal_data, year, watershed) {
  
  # Convert DOY to date for better x-axis labels
  doy_to_date <- function(doy, year = 2024) {
    as.Date(doy - 1, origin = paste0(year, "-01-01"))
  }
  
  # Fixed x-axis range: June 1 (DOY 152) to July 30 (DOY 211)
  june_1_doy <- 152  # June 1
  july_30_doy <- 211  # July 30
  
  # Create custom x-axis breaks every 10 days
  doy_breaks <- seq(june_1_doy, july_30_doy, by = 10)
  
  ggplot(natal_data, aes(x = DOY, y = dailyCPUEprop)) +
    geom_line(color = "steelblue", linewidth = 1.5, alpha = 0.8) +
    geom_ribbon(aes(ymin = 0, ymax = dailyCPUEprop), 
                fill = "steelblue", alpha = 0.3) +
    scale_x_continuous(
      breaks = doy_breaks,
      labels = function(x) {
        paste0(x, "\n", format(doy_to_date(x), "%b %d"))
      },
      limits = c(june_1_doy, july_30_doy),  # Fixed limits for all years
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(x * 100, 1), "%"),
      limits = c(0, max(natal_data$dailyCPUEprop, na.rm = TRUE) * 1.05)
    ) +
    labs(
      title = paste("Daily CPUE Distribution -", year),
      subtitle = paste(watershed, "Watershed | Fish catch timing throughout the season"),
      x = "Day of Year (Date)",
      y = "Daily CPUE Proportion"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
      axis.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

#' Create enhanced cumulative distribution plots with FIXED management unit names
create_enhanced_cumulative_plots <- function(cumulative_data, summary_stats, 
                                             watershed, interval_days) {
  
  cat("Creating enhanced visualization plots with FIXED management unit names...\n")
  
  # Create output directory
  plot_dir <- here("Analysis_Results/Cumulative_Distribution/Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # CORRECTED MANAGEMENT UNIT ORDER - Based on actual data from your scripts
  # These are the ACTUAL names found in your data
  desired_order <- c(
    "N. Fork Kusko",
    "E. Fork Kuskokwim",        # No "River" - this appears to be the actual name
    "S. Fork Kusko", 
    "Upper Kusko Main",
    "Big River",
    "Takotna and Nixon Fork",
    "Tatlawiksuk",
    "Swift",
    "Stony", 
    "Holitna and Hoholitna",    # This appears to be the actual combined name
    "Middle Kusko Main",
    "George",
    "Oskakawlik", 
    "Holokuk",
    "Aniak",
    "Tuluksak",
    "Kisaralik",
    "Kwethluk",
    "Johnson",
    "Lower Kusko"
  )
  
  # Get actual management units from data
  actual_mgmt_rivers <- sort(unique(cumulative_data$mgmt_river))
  cat(paste("ACTUAL management units in data:", length(actual_mgmt_rivers), "\n"))
  for (i in 1:length(actual_mgmt_rivers)) {
    cat(paste("  ", i, ". '", actual_mgmt_rivers[i], "'\n", sep = ""))
  }
  
  # Use ONLY the management units that actually exist in the data
  # Don't try to impose our desired order if the names don't match
  mgmt_rivers_for_plotting <- actual_mgmt_rivers
  
  # Try to match as many as possible from desired order, but include all actual units
  matched_units <- intersect(desired_order, actual_mgmt_rivers)
  unmatched_units <- setdiff(actual_mgmt_rivers, desired_order)
  
  cat(paste("MATCHED units from desired order:", length(matched_units), "\n"))
  cat(paste("UNMATCHED units (will be added at end):", length(unmatched_units), "\n"))
  for (unit in unmatched_units) {
    cat(paste("  UNMATCHED: '", unit, "'\n", sep = ""))
  }
  
  # Create final ordering: matched units in desired order + unmatched units at end
  final_order <- c(matched_units, unmatched_units)
  
  cat("\nFINAL ORDER FOR PLOTTING:\n")
  for (i in 1:length(final_order)) {
    cat(paste("  ", i, ". '", final_order[i], "'\n", sep = ""))
  }
  
  # Create red → orange → blue gradient function
  create_watershed_colors <- function(mgmt_rivers_ordered) {
    n_watersheds <- length(mgmt_rivers_ordered)
    
    # Create gradient from red → orange → blue
    color_gradient <- colorRampPalette(c(
      "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
      "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
    ))(n_watersheds)
    
    # Create named vector
    names(color_gradient) <- mgmt_rivers_ordered
    return(color_gradient)
  }
  
  # Generate colors for actual management units
  watershed_colors <- create_watershed_colors(final_order)
  
  # Create ordered factor for proper legend ordering
  cumulative_data$mgmt_river <- factor(cumulative_data$mgmt_river, levels = final_order)
  
  cat(paste("Applied colors to", length(watershed_colors), "management units\n"))
  
  #--------------------------------------------------------------------------#
  # CREATE INDIVIDUAL PLOTS FOR EACH YEAR
  #--------------------------------------------------------------------------#
  
  individual_plots <- list()
  years <- sort(unique(cumulative_data$year))
  
  for (year in years) {
    year_data <- cumulative_data %>% filter(year == !!year)
    
    # Load natal data for CPUE histogram
    natal_data <- load_natal_data(year, watershed)
    
    # Create CPUE histogram for this year
    cpue_histogram <- create_cpue_histogram(natal_data, year, watershed)
    
    # Get ALL unique DOY values for this year
    all_doys_this_year <- sort(unique(year_data$doy))
    
    # Create enhanced plot with NO POINTS, clean lines only
    p_year <- ggplot(year_data, aes(x = doy, y = cumulative_percent_of_unit, color = mgmt_river)) +
      # Add reference lines for key timing milestones
      geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.7) +
      geom_hline(yintercept = c(10, 90), color = "gray80", linetype = "dotted", alpha = 0.5) +
      
      # Main data - LINES ONLY (NO POINTS)
      geom_line(linewidth = 2.5, alpha = 0.9) +
      
      # Color and scales
      scale_color_manual(values = watershed_colors, 
                         breaks = final_order,
                         guide = guide_legend(override.aes = list(linewidth = 3))) +
      scale_y_continuous(
        labels = function(x) paste0(round(x, 1), "%"),
        limits = c(0, 100),
        breaks = c(0, 10, 25, 50, 75, 90, 100),
        minor_breaks = seq(0, 100, 5)
      ) +
      # Show all DOY values in the data
      scale_x_continuous(
        breaks = all_doys_this_year,
        labels = function(x) {
          date_str <- format(as.Date(x - 1, origin = paste0(year, "-01-01")), "%b %d")
          paste0("DOY ", x, "\n", date_str)
        },
        limits = c(min(year_data$doy), max(year_data$doy)),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      
      # Enhanced labels and theme
      labs(
        title = paste("FIXED: Salmon Run Timing Progress -", year),
        subtitle = paste(watershed, "Watershed | FIXED management unit names | ALL", length(all_doys_this_year), "dates shown"),
        x = paste("Day of Year (Date) - All", interval_days, "Day Intervals Shown"),
        y = "Cumulative Percent of Management Unit's Total Production",
        color = "Management Unit\n(Watershed Position)",
        caption = paste("FIXED: Corrected management unit names | No points, clean lines |", length(all_doys_this_year), "data points | Red=upstream, Blue=downstream")
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "gray20"),
        plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 10, hjust = 0.5, color = "gray50"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        axis.line = element_line(color = "gray20", linewidth = 1.2),
        axis.text = element_text(size = 10, color = "gray20"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, color = "gray20"),
        axis.title = element_text(face = "bold", size = 12, color = "gray20"),
        axis.ticks = element_line(color = "gray20", linewidth = 0.8),
        panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
        panel.border = element_rect(color = "gray20", fill = NA, linewidth = 1.0),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    # Save enhanced combined plot
    combined_filename <- paste0("FIXED_run_timing_progress_", year, "_", watershed, "_corrected_names.png")
    
    png(file.path(plot_dir, combined_filename), width = 18, height = 12, units = "in", res = 300, bg = "white")
    grid.arrange(p_year, cpue_histogram, nrow = 2, heights = c(3, 1))
    dev.off()
    
    individual_plots[[as.character(year)]] <- p_year
    cat(paste("✓ Created FIXED timing plot for", year, "with corrected names | File:", combined_filename, "\n"))
  }
  
  #--------------------------------------------------------------------------#
  # CREATE OVERVIEW FACETED PLOT - ALL YEARS TOGETHER
  #--------------------------------------------------------------------------#
  
  # Order data for consistent legend ordering
  cumulative_data$mgmt_river <- factor(cumulative_data$mgmt_river, levels = final_order)
  
  # Use all DOY values present in the data for the overview plot
  all_doys_in_overview <- sort(unique(cumulative_data$doy))
  doy_range_overview <- range(cumulative_data$doy)
  
  cat(paste("Overview plot: showing", length(all_doys_in_overview), "unique dates across all years\n"))
  
  # Paper-optimized faceted plot with enhanced readability
  p1_faceted <- ggplot(cumulative_data, aes(x = doy, y = cumulative_percent_of_unit, color = mgmt_river)) +
    geom_line(linewidth = 2.0, alpha = 0.9) +  # Thicker lines for paper
    facet_wrap(~year, ncol = 1, scales = "free_x") +
    scale_color_manual(values = watershed_colors,
                       breaks = final_order,
                       guide = guide_legend(override.aes = list(linewidth = 3), # Thicker legend lines
                                            ncol = 1,
                                            keywidth = 1.5,
                                            keyheight = 1.0)) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, 100),
      breaks = seq(0, 100, 25)
    ) +
    scale_x_continuous(
      limits = doy_range_overview,
      breaks = seq(from = min(all_doys_in_overview), to = max(all_doys_in_overview), by = 9),
      labels = function(x) {
        date_obj <- as.Date(x - 1, origin = "2024-01-01")
        format(date_obj, "%b %d")
      },
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = paste("Cumulative Distribution of returning Chinook stocks over time"),
      x = "Date",
      y = "Cumulative Percent of Unit's Total Production",
      color = "Management Unit"
    ) +
    theme_minimal(base_size = 14) +  # Larger base font for paper
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 14),
      
      # Enhanced legend for paper
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 11),
      legend.margin = margin(l = 20),
      legend.box.spacing = unit(0.5, "cm"),
      
      # Enhanced axis styling for paper readability
      axis.line = element_line(color = "gray20", linewidth = 1.0),
      axis.text = element_text(color = "gray20", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "gray20"),
      axis.title = element_text(face = "bold", color = "gray20", size = 13),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      
      # Panel optimizations
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.4),
      panel.spacing = unit(1.5, "lines"),
      
      # Generous margins for paper
      plot.margin = margin(15, 15, 15, 15, "mm")
    )
  
  # Save with paper-optimized dimensions
  ggsave(file.path(plot_dir, paste0("FIXED_cumulative_overview_", watershed, "_corrected_names_PAPER.png")), 
         p1_faceted, 
         width = 16, height = 20, 
         dpi = 300, bg = "white")
  
  cat("✓ Created FIXED overview faceted plot with corrected names\n")
  
  cat("FIXED plots created with corrected management unit names!\n")
  
  return(list(
    individual_year_plots = individual_plots,
    overview_faceted = p1_faceted,
    mgmt_units_used = final_order,
    colors_used = watershed_colors
  ))
}

#' Create horizontal range plot showing min-max with mean marked
create_run_duration_range_plot <- function(run_duration_data, watershed_colors, mgmt_rivers_ordered, 
                                           watershed, interval_days) {
  
  cat("Creating horizontal range plot with mean markers...\n")
  
  # Create output directory
  plot_dir <- here("Analysis_Results/Cumulative_Distribution/Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Clean data: only keep exactly the management units that are in mgmt_rivers_ordered
  clean_data <- run_duration_data %>%
    filter(!is.na(mgmt_river), !is.na(days_to_50_percent)) %>%
    filter(mgmt_river %in% mgmt_rivers_ordered)
  
  cat(paste("Data filtering:\n"))
  cat(paste("  Original run duration rows:", nrow(run_duration_data), "\n"))
  cat(paste("  After cleaning:", nrow(clean_data), "\n"))
  cat(paste("  Management units:", length(unique(clean_data$mgmt_river)), "\n"))
  
  if (nrow(clean_data) == 0) {
    cat("❌ No clean data available for range plot\n")
    return(list(range_plot = NULL, summary_stats = NULL, filename = NULL))
  }
  
  # Ensure clean_data has the same factor levels as the cumulative plots
  clean_data$mgmt_river <- factor(clean_data$mgmt_river, levels = mgmt_rivers_ordered)
  
  # Calculate summary statistics for the range plot
  duration_stats <- clean_data %>%
    group_by(mgmt_river) %>%
    summarise(
      n_years = n(),
      mean_days = mean(days_to_50_percent, na.rm = TRUE),
      min_days = min(days_to_50_percent, na.rm = TRUE),
      max_days = max(days_to_50_percent, na.rm = TRUE),
      range_days = max_days - min_days,
      .groups = 'drop'
    ) %>%
    arrange(mean_days)
  
  # Ensure factor levels match for plotting
  duration_stats$mgmt_river <- factor(duration_stats$mgmt_river, levels = mgmt_rivers_ordered)
  
  # Create the MODERN HORIZONTAL range plot
  p_range <- ggplot(duration_stats, aes(y = mgmt_river)) +
    # Range lines (min to max)
    geom_segment(aes(x = min_days, xend = max_days, 
                     color = mgmt_river), 
                 linewidth = 6, alpha = 0.7) +
    
    # Mean markers (black marks)
    geom_point(aes(x = mean_days), 
               color = "black", 
               size = 4, 
               shape = "|",  # Vertical line shape
               stroke = 2) +
    
    # Use the exact same colors as cumulative plots for the range lines
    scale_color_manual(values = watershed_colors, 
                       breaks = mgmt_rivers_ordered,
                       guide = "none") +
    
    # Modern X-axis formatting
    scale_x_continuous(
      name = "Days to 50% return",
      breaks = function(x) pretty(x, n = 6),
      expand = expansion(mult = c(0.02, 0.02)),
      labels = function(x) paste0(x, " days")
    ) +
    
    # Clean Y-axis formatting
    scale_y_discrete(
      name = NULL,  # Remove y-axis title for cleaner look
      limits = rev(mgmt_rivers_ordered)
    ) +
    
    # Modern, minimal labels
    labs(
      title = "Days to 50% of return by stock"
    ) +
    
    # Modern minimal theme
    theme_void() +
    theme(
      # Modern typography
      plot.title = element_text(
        face = "bold", 
        size = 20, 
        hjust = 0, 
        color = "gray15",
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = 13, 
        hjust = 0, 
        color = "gray50",
        margin = margin(b = 20)
      ),
      plot.caption = element_text(
        size = 11, 
        hjust = 0, 
        color = "gray60",
        margin = margin(t = 15)
      ),
      
      # Clean axis styling
      axis.text.y = element_text(
        size = 11, 
        color = "gray30",
        hjust = 1,
        margin = margin(r = 10)
      ),
      axis.text.x = element_text(
        size = 10, 
        color = "gray30",
        margin = margin(t = 8)
      ),
      axis.title.x = element_text(
        size = 12, 
        color = "gray20",
        face = "bold",
        margin = margin(t = 15)
      ),
      
      # Minimal grid
      panel.grid.major.x = element_line(
        color = "gray92", 
        linewidth = 0.4,
        linetype = "solid"
      ),
      panel.grid.minor.x = element_line(
        color = "gray96", 
        linewidth = 0.2
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      # Modern spacing
      plot.margin = margin(20, 25, 20, 20, "mm"),
      
      # Clean background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Remove all axis lines and ticks for ultra-clean look
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Save the modern range plot
  range_filename <- paste0("FIXED_run_duration_range_", watershed, "_", interval_days, "day_intervals.png")
  ggsave(file.path(plot_dir, range_filename), p_range, 
         width = 12, height = 8, dpi = 300, bg = "white")
  
  cat(paste("✅ Created FIXED modern horizontal range plot:", range_filename, "\n"))
  cat(paste("   Showing min-max range with mean for", length(unique(clean_data$mgmt_river)), "management units\n"))
  cat(paste("   Colors match cumulative distribution plots EXACTLY\n"))
  
  # Print summary statistics
  cat("\nRun Duration Range Statistics:\n")
  print(duration_stats, n = Inf)
  
  return(list(
    range_plot = p_range,
    summary_stats = duration_stats,
    filename = range_filename
  ))
}

#' Export enhanced cumulative data
export_enhanced_cumulative_data <- function(cumulative_data, summary_stats, run_duration_data, 
                                            watershed, interval_days) {
  
  cat("Exporting data with corrected management unit names...\n")
  
  # Create output directory
  csv_dir <- here("Analysis_Results/Cumulative_Distribution/CSV")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Export full cumulative data
  cumulative_clean <- cumulative_data %>%
    select(year, doy, interval_number, mgmt_river, cumulative_production, 
           cumulative_percent_of_unit, overall_percent, fish_count_cumulative) %>%
    arrange(year, mgmt_river, doy)
  
  cumulative_filepath <- file.path(csv_dir, paste0("FIXED_cumulative_distribution_data_", 
                                                   watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(cumulative_clean, cumulative_filepath, row.names = FALSE)
  cat(paste("Exported FIXED cumulative data to:", basename(cumulative_filepath), "\n"))
  
  # Export run duration metrics
  run_duration_clean <- run_duration_data %>%
    select(year, mgmt_river, run_start_doy, doy_above_50_percent, days_to_above_50_percent, 
           actual_percent_at_doy, total_production, final_overall_percent) %>%
    arrange(year, mgmt_river)
  
  run_duration_filepath <- file.path(csv_dir, paste0("FIXED_run_duration_metrics_", 
                                                     watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(run_duration_clean, run_duration_filepath, row.names = FALSE)
  cat(paste("Exported FIXED run duration metrics to:", basename(run_duration_filepath), "\n"))
  
  # Export summary statistics
  summary_filepath <- file.path(csv_dir, paste0("FIXED_timing_summary_stats_", 
                                                watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(summary_stats, summary_filepath, row.names = FALSE)
  cat(paste("Exported FIXED summary statistics to:", basename(summary_filepath), "\n"))
  
  # Export wide format for easy plotting
  cumulative_wide <- cumulative_data %>%
    select(year, doy, mgmt_river, cumulative_percent_of_unit) %>%
    pivot_wider(names_from = mgmt_river, values_from = cumulative_percent_of_unit, 
                values_fill = 0, names_prefix = "mgmt_") %>%
    arrange(year, doy)
  
  wide_filepath <- file.path(csv_dir, paste0("FIXED_cumulative_distribution_wide_", 
                                             watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(cumulative_wide, wide_filepath, row.names = FALSE)
  cat(paste("Exported FIXED wide format data to:", basename(wide_filepath), "\n"))
  
  cat("All FIXED CSV exports completed!\n")
}

#' Calculate cumulative distributions by management unit with 3-day intervals
analyze_cumulative_distributions <- function(years = c(2017, 2018, 2019, 2020, 2021),
                                             watershed = "Kusko",
                                             interval_days = 3,
                                             export_csv = TRUE) {
  
  cat("=== Starting FIXED Cumulative Distribution Analysis ===\n")
  cat(paste("Watershed:", watershed, "\n"))
  cat(paste("Years:", paste(years, collapse = ", "), "\n"))
  cat(paste("Data points every", interval_days, "days\n"))
  cat("FIXED: Management unit names corrected to match actual data\n")
  
  # Get watershed parameters
  if (watershed == "Kusko") {
    params <- list(sensitivity_threshold = 0.7, min_error = 0.0006, min_stream_order = 3)
  } else if (watershed == "Yukon") {
    params <- list(sensitivity_threshold = 0.7, min_error = 0.003, min_stream_order = 5)
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon'")
  }
  
  # Storage for all results
  all_cumulative_data <- list()
  all_summary_stats <- list()
  all_run_duration_data <- list()
  
  ################################################################################
  # PROCESS EACH YEAR
  ################################################################################
  
  for (year in years) {
    cat(paste("\n--- Processing year", year, "---\n"))
    
    # Load spatial data
    spatial_data <- load_spatial_data(watershed, 8, params$min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    
    # Load natal data
    natal_data <- load_natal_data(year, watershed)
    cat(paste("Loaded", nrow(natal_data), "fish observations\n"))
    
    # Check for management units
    if (!"mgmt_river" %in% colnames(edges)) {
      warning("Management river data not found. Skipping year ", year)
      next
    }
    
    # DIAGNOSTIC: Print actual management unit names found in data
    actual_mgmt_units <- sort(unique(edges$mgmt_river[!is.na(edges$mgmt_river) & edges$mgmt_river != ""]))
    cat(paste("FOUND", length(actual_mgmt_units), "management units in spatial data for", year, ":\n"))
    for (i in 1:length(actual_mgmt_units)) {
      cat(paste("  ", i, ". '", actual_mgmt_units[i], "'\n", sep = ""))
    }
    
    # Setup assignment parameters
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, params$min_error)
    priors <- setup_watershed_priors(edges, params$min_stream_order, watershed, natal_data)
    
    # Determine DOY range and create 3-day intervals
    doy_range <- range(natal_data$DOY, na.rm = TRUE)
    min_doy <- floor(doy_range[1])
    max_doy <- ceiling(doy_range[2])
    
    # Create sequence of DOY breakpoints every 3 days starting from min_doy
    doy_breakpoints <- seq(from = min_doy, to = max_doy, by = interval_days)
    
    # Ensure we include the final DOY if it's not exactly on a 3-day boundary
    if (max(doy_breakpoints) < max_doy) {
      doy_breakpoints <- c(doy_breakpoints, max_doy)
    }
    
    cat(paste("DOY range:", min_doy, "to", max_doy, "\n"))
    cat(paste("Created", length(doy_breakpoints), "breakpoints every", interval_days, "days\n"))
    
    # Calculate cumulative distributions at each 3-day interval
    year_cumulative_data <- data.frame()
    
    for (i in 1:length(doy_breakpoints)) {
      current_doy <- doy_breakpoints[i]
      
      # Filter data up to current DOY (cumulative)
      cumulative_data_subset <- natal_data %>%
        filter(DOY <= current_doy)
      
      if (nrow(cumulative_data_subset) == 0) {
        next
      }
      
      cat(paste("Processing DOY", current_doy, "with", nrow(cumulative_data_subset), "fish\n"))
      
      # Perform assignment for cumulative data up to this DOY
      assignment_matrix <- perform_assignment(
        cumulative_data_subset, edges, watershed, priors, pid_iso, error, params$sensitivity_threshold
      )
      
      # Calculate stream-level production
      stream_production <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
      
      # Aggregate by management units
      mgmt_result <- process_mgmt_river_data(edges, stream_production)
      
      if (!is.null(mgmt_result)) {
        # Remove spatial geometry if present
        if (inherits(mgmt_result, "sf")) {
          mgmt_result <- st_drop_geometry(mgmt_result)
        }
        
        # Calculate total production across all units
        total_production_all <- sum(mgmt_result$total_production, na.rm = TRUE)
        
        # Create data frame for this DOY interval
        doy_data <- mgmt_result %>%
          mutate(
            year = year,
            doy = current_doy,
            interval_number = i,
            cumulative_production = total_production,
            fish_count_cumulative = nrow(cumulative_data_subset),
            total_production_all_units = total_production_all
          ) %>%
          select(year, doy, interval_number, mgmt_river, cumulative_production, 
                 fish_count_cumulative, total_production_all_units, edge_count)
        
        year_cumulative_data <- bind_rows(year_cumulative_data, doy_data)
      }
    }
    
    if (nrow(year_cumulative_data) > 0) {
      # For each management unit, calculate final production (at end of season)
      final_production_by_unit <- year_cumulative_data %>%
        group_by(mgmt_river) %>%
        reframe(final_production = max(cumulative_production, na.rm = TRUE))
      
      # Now calculate cumulative percentages relative to each unit's final production
      year_cumulative_data <- year_cumulative_data %>%
        left_join(final_production_by_unit, by = "mgmt_river") %>%
        mutate(
          # Key metric: What % of THIS UNIT'S final production has been reached?
          cumulative_percent_of_unit = (cumulative_production / final_production) * 100,
          # Also keep overall proportion for reference
          overall_proportion = cumulative_production / total_production_all_units,
          overall_percent = (cumulative_production / total_production_all_units) * 100
        ) %>%
        select(-final_production)  # Remove temporary column
      
      # Calculate summary statistics for this year
      year_summary <- year_cumulative_data %>%
        group_by(mgmt_river) %>%
        arrange(doy) %>%
        mutate(
          final_production = max(cumulative_production, na.rm = TRUE),
          final_overall_percent = max(overall_percent, na.rm = TRUE)
        ) %>%
        reframe(
          year = first(year),
          total_production = max(cumulative_production, na.rm = TRUE),
          final_overall_percent = max(overall_percent, na.rm = TRUE),
          
          # Timing statistics (DOY when unit reaches certain percentages of its OWN final production)
          doy_10_percent = doy[which.min(abs(cumulative_percent_of_unit - 10))],
          doy_25_percent = doy[which.min(abs(cumulative_percent_of_unit - 25))],
          doy_50_percent = doy[which.min(abs(cumulative_percent_of_unit - 50))],
          doy_75_percent = doy[which.min(abs(cumulative_percent_of_unit - 75))],
          doy_90_percent = doy[which.min(abs(cumulative_percent_of_unit - 90))],
          
          # Timing spread
          timing_spread = doy_90_percent - doy_10_percent,
          median_timing = doy_50_percent
        )
      
      # Calculate run duration metrics - SIMPLE: First date above 50%
      year_run_duration <- year_cumulative_data %>%
        group_by(mgmt_river) %>%
        arrange(doy) %>%
        reframe(
          year = first(year),
          run_start_doy = min(doy, na.rm = TRUE),  # First day with any fish
          
          # SIMPLE: First DOY where cumulative percent goes above 50%
          doy_above_50_percent = {
            above_50_idx <- which(cumulative_percent_of_unit > 50)
            if (length(above_50_idx) > 0) {
              doy[min(above_50_idx)]  # First DOY above 50%
            } else {
              max(doy, na.rm = TRUE)  # If never reaches 50%, use last DOY
            }
          },
          
          # Actual percentage at that DOY
          actual_percent_at_doy = {
            above_50_idx <- which(cumulative_percent_of_unit > 50)
            if (length(above_50_idx) > 0) {
              cumulative_percent_of_unit[min(above_50_idx)]
            } else {
              max(cumulative_percent_of_unit, na.rm = TRUE)
            }
          },
          
          # Days from run start to above 50%
          days_to_above_50_percent = doy_above_50_percent - run_start_doy,
          
          # Additional info
          max_percent_reached = max(cumulative_percent_of_unit, na.rm = TRUE),
          total_production = max(cumulative_production, na.rm = TRUE),
          final_overall_percent = max(overall_percent, na.rm = TRUE)
        )
      
      # Print simple debugging info for first year
      if (year == years[1]) {
        cat("\n=== SIMPLE 50% CALCULATION FOR", year, "===\n")
        debug_sample <- year_run_duration %>% 
          select(mgmt_river, run_start_doy, doy_above_50_percent, actual_percent_at_doy, 
                 days_to_above_50_percent, max_percent_reached) %>%
          head(5)
        print(debug_sample)
      }
      
      all_run_duration_data[[as.character(year)]] <- year_run_duration
      all_summary_stats[[as.character(year)]] <- year_summary
      all_cumulative_data[[as.character(year)]] <- year_cumulative_data
    }
    
    cat(paste("Completed year", year, "with", nrow(year_cumulative_data), "data points\n"))
  }
  
  ################################################################################
  # COMBINE ALL YEARS AND CREATE PLOTS
  ################################################################################
  
  if (length(all_cumulative_data) == 0) {
    stop("No valid data processed for any year")
  }
  
  # Combine all data
  combined_cumulative_data <- bind_rows(all_cumulative_data)
  combined_summary_stats <- bind_rows(all_summary_stats)
  combined_run_duration <- bind_rows(all_run_duration_data)
  
  cat(paste("Combined data:", nrow(combined_cumulative_data), "total data points\n"))
  cat(paste("Management units analyzed:", length(unique(combined_cumulative_data$mgmt_river)), "\n"))
  
  # Print actual management units found in final data
  final_mgmt_units <- sort(unique(combined_cumulative_data$mgmt_river))
  cat("\nFINAL MANAGEMENT UNITS IN COMBINED DATA:\n")
  for (i in 1:length(final_mgmt_units)) {
    cat(paste("  ", i, ". '", final_mgmt_units[i], "'\n", sep = ""))
  }
  
  # Create enhanced visualization plots
  plots <- create_enhanced_cumulative_plots(combined_cumulative_data, combined_summary_stats, 
                                            watershed, interval_days)
  
  # CREATE RUN DURATION BOXPLOT WITH SAME MANAGEMENT UNITS AND COLORS
  cat("\n=== Creating Run Duration Boxplot with FIXED Management Units ===\n")
  
  # Use the EXACT same management units and colors as the cumulative plots
  final_order_from_plots <- plots$mgmt_units_used
  watershed_colors_from_plots <- plots$colors_used
  
  cat(paste("Using", length(final_order_from_plots), "management units for boxplot\n"))
  
  # Filter run duration data to match the cumulative plots exactly
  matching_run_duration <- combined_run_duration %>%
    filter(mgmt_river %in% final_order_from_plots) %>%
    filter(!is.na(days_to_above_50_percent))  # Use the simple above 50% calculation
  
  if (nrow(matching_run_duration) > 0) {
    # Use simple above 50% values for the boxplot
    matching_run_duration$days_to_50_percent <- matching_run_duration$days_to_above_50_percent
    
    boxplot_result <- create_run_duration_range_plot(matching_run_duration, watershed_colors_from_plots, 
                                                     final_order_from_plots, watershed, interval_days)
    
    # Add range plot to plots list
    plots$run_duration_range_plot <- boxplot_result$range_plot
    plots$range_plot_summary_stats <- boxplot_result$summary_stats
    
    cat("✓ Successfully created run duration RANGE PLOT with SIMPLE above 50% calculations\n")
  } else {
    cat("❌ No run duration data available for range plot\n")
    plots$run_duration_range_plot <- NULL
    plots$range_plot_summary_stats <- NULL
  }
  
  # Export data if requested
  if (export_csv) {
    export_enhanced_cumulative_data(combined_cumulative_data, combined_summary_stats, 
                                    combined_run_duration, watershed, interval_days)
  }
  
  return(list(
    cumulative_data = combined_cumulative_data,
    summary_stats = combined_summary_stats,
    run_duration_metrics = combined_run_duration,
    plots = plots,
    parameters = list(
      watershed = watershed,
      years = years,
      interval_days = interval_days,
      params = params
    )
  ))
}

################################################################################
# MAIN EXECUTION - RUNS AUTOMATICALLY WHEN SCRIPT IS EXECUTED
################################################################################

cat("################################################################################\n")
cat("# EXECUTING FIXED CUMULATIVE DISTRIBUTION ANALYSIS\n")
cat("################################################################################\n")
cat("FIXES APPLIED:\n")
cat("  - Diagnostic printing of actual management unit names found in data\n")
cat("  - Use actual names instead of imposing desired order with wrong names\n")
cat("  - Match as many units as possible from desired order, add unmatched at end\n")
cat("  - All plots and exports will use 'FIXED_' prefix to distinguish from old versions\n")
cat("################################################################################\n")

# Run the FIXED analysis
cat("Starting FIXED Cumulative Distribution Analysis...\n")
results <- tryCatch({
  analyze_cumulative_distributions(
    years = c(2017, 2018, 2019, 2020, 2021),
    watershed = "Kusko",
    interval_days = 3,
    export_csv = TRUE
  )
}, error = function(e) {
  cat("Error during analysis:", e$message, "\n")
  cat("Please check that all required files and paths are correct.\n")
  return(NULL)
})

if (!is.null(results)) {
  cat("\n=== FIXED ANALYSIS COMPLETE ===\n")
  cat(paste("Total data points:", nrow(results$cumulative_data), "\n"))
  cat(paste("Management units analyzed:", length(unique(results$cumulative_data$mgmt_river)), "\n"))
  cat(paste("DOY range:", min(results$cumulative_data$doy), "to", max(results$cumulative_data$doy), "\n"))
  cat(paste("Data points created every", results$parameters$interval_days, "days\n"))
  
  # Show the final management units that were actually used
  final_mgmt_units <- unique(results$cumulative_data$mgmt_river)
  cat("\n=== FINAL MANAGEMENT UNITS USED (IN ORDER) ===\n")
  for (i in 1:length(final_mgmt_units)) {
    cat(paste("  ", i, ". '", final_mgmt_units[i], "'\n", sep = ""))
  }
  
  # Show run duration metrics if available
  if (!is.null(results$run_duration_metrics) && nrow(results$run_duration_metrics) > 0) {
    cat("\n=== RUN DURATION METRICS ===\n")
    cat(paste("Run duration data calculated for", nrow(results$run_duration_metrics), "stock-year combinations\n"))
    
    # Calculate summary statistics - SIMPLE
    duration_summary <- results$run_duration_metrics %>%
      group_by(mgmt_river) %>%
      reframe(
        avg_days_to_above_50 = round(mean(days_to_above_50_percent, na.rm = TRUE), 1),
        min_days = min(days_to_above_50_percent, na.rm = TRUE),
        max_days = max(days_to_above_50_percent, na.rm = TRUE),
        avg_actual_percent = round(mean(actual_percent_at_doy, na.rm = TRUE), 1)
      ) %>%
      arrange(avg_days_to_above_50)
    
    cat("\nAverage days from run start to ABOVE 50% completion by management unit:\n")
    print(duration_summary)
    
    # Show which management unit has fastest/slowest average build-up
    if (nrow(duration_summary) > 0) {
      fastest_unit <- duration_summary$mgmt_river[1]
      slowest_unit <- duration_summary$mgmt_river[nrow(duration_summary)]
      fastest_days <- duration_summary$avg_days_to_above_50[1]
      slowest_days <- duration_summary$avg_days_to_above_50[nrow(duration_summary)]
      
      cat(paste("\n🏃 FASTEST BUILD-UP:", fastest_unit, "- avg", fastest_days, "days to above 50%\n"))
      cat(paste("🐌 SLOWEST BUILD-UP:", slowest_unit, "- avg", slowest_days, "days to above 50%\n"))
    }
  }
  
  # Verify 3-day intervals
  doy_intervals <- results$cumulative_data %>%
    filter(mgmt_river == first(mgmt_river)) %>%
    filter(year == first(year)) %>%
    arrange(doy) %>%
    mutate(interval_gap = doy - lag(doy, default = first(doy))) %>%
    filter(row_number() > 1)
  
  cat("\n=== INTERVAL VERIFICATION ===\n")
  cat("Gaps between consecutive data points (should be 3 days):\n")
  print(table(doy_intervals$interval_gap))
  
  # Show files created
  years_processed <- sort(unique(results$cumulative_data$year))
  cat("\n=== FIXED FILES CREATED ===\n")
  cat("📊 FIXED INDIVIDUAL YEAR TIMING PLOTS:\n")
  for (year in years_processed) {
    filename <- paste0("FIXED_run_timing_progress_", year, "_", results$parameters$watershed, "_corrected_names.png")
    cat(paste("   ", year, ":", filename, "\n"))
  }
  
  cat("\n📊 FIXED OVERVIEW PLOT:\n")
  cat(paste("   FIXED_cumulative_overview_", results$parameters$watershed, "_corrected_names_PAPER.png\n"))
  
  cat("\n📊 FIXED RUN DURATION RANGE PLOT:\n")
  if (!is.null(results$plots$run_duration_range_plot)) {
    cat(paste("   FIXED_run_duration_range_", results$parameters$watershed, "_", results$parameters$interval_days, "day_intervals.png\n"))
    cat("   → ✨ MODERN DESIGN: Clean range plot with mean markers\n")
    cat("   → FIXED: Horizontal orientation (left to right)\n")
    cat("   → FIXED: Uses corrected management unit names\n")
    cat("   → FIXED: Colors match exactly with cumulative distribution plots\n")
    cat("   → NEW: Shows min-max range with black mark at mean (no boxplot)\n")
  } else {
    cat("   ❌ Run duration range plot not created (no valid data)\n")
  }
  
  cat("\n📁 FIXED CSV FILES:\n")
  cat("   ✓ FIXED_cumulative_distribution_data_[watershed]_[interval]day_intervals.csv\n")
  cat("   ✓ FIXED_run_duration_metrics_[watershed]_[interval]day_intervals.csv\n") 
  cat("   ✓ FIXED_timing_summary_stats_[watershed]_[interval]day_intervals.csv\n")
  cat("   ✓ FIXED_cumulative_distribution_wide_[watershed]_[interval]day_intervals.csv\n")
  
  cat(paste("\n📁 All files location:", here("Analysis_Results/Cumulative_Distribution/"), "\n"))
  
  cat("\n🎯 FIXES APPLIED:\n")
  cat("   ✅ DIAGNOSTIC: Prints actual management unit names found in data\n")
  cat("   ✅ FLEXIBLE ORDERING: Uses actual names, doesn't force incorrect desired names\n")
  cat("   ✅ SMART MATCHING: Matches what it can from desired order, includes all actual units\n")
  cat("   ✅ CLEAR NAMING: All output files have 'FIXED_' prefix\n")
  cat("   ✅ NO MORE MISSING RIVERS: All management units in data will be included\n")
  
  cat("\n✅ FIXED: All management rivers now showing up correctly!\n")
  cat("Fixed analysis completed successfully!\n")
  
} else {
  cat("Analysis failed. Please check error messages above.\n")
}

cat("################################################################################\n")
cat("# FIXED SCRIPT EXECUTION COMPLETE\n")
cat("################################################################################\n")

