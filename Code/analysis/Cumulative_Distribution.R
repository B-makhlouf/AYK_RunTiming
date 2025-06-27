
################################################################################
# COMPLETE CUMULATIVE DISTRIBUTION ANALYSIS SCRIPT
################################################################################
# PURPOSE: Calculate cumulative distribution plots for each management unit
#          with data points every 3 days - READY TO RUN
# USAGE:   Simply execute this script - no sourcing required
# OUTPUT:  Plots and CSV data showing how fish production accumulates over time
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
})

# Source required utilities
cat("Loading spatial utilities and assignment functions...\n")
tryCatch({
  source(here("code/utils/spatial_utils.R"))
  source(here("code/assignment.R"))
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

#' Combine timing plot with CPUE histogram
create_combined_timing_plot <- function(timing_plot, cpue_histogram, year, watershed) {
  
  # Create a function that draws both plots
  combined_plot_function <- function() {
    # Set up layout with grid
    grid.newpage()
    
    # Create layout: 70% for timing plot, 30% for histogram
    pushViewport(viewport(layout = grid.layout(3, 1, 
                                               heights = unit(c(0.1, 0.6, 0.3), "npc"))))
    
    # Add overall title
    pushViewport(viewport(layout.pos.row = 1))
    grid.text(paste("Salmon Run Analysis -", year, "|", watershed, "Watershed"),
              gp = gpar(fontsize = 16, fontface = "bold", col = "gray20"))
    popViewport()
    
    # Add timing plot
    pushViewport(viewport(layout.pos.row = 2))
    print(timing_plot, newpage = FALSE)
    popViewport()
    
    # Add CPUE histogram
    pushViewport(viewport(layout.pos.row = 3))
    print(cpue_histogram, newpage = FALSE)
    popViewport()
    
    popViewport()
  }
  
  return(combined_plot_function)
}

#' Calculate cumulative distributions by management unit with 3-day intervals
analyze_cumulative_distributions <- function(years = c(2017, 2018, 2019, 2020, 2021),
                                             watershed = "Kusko",
                                             interval_days = 3,
                                             export_csv = TRUE) {
  
  cat("=== Starting Cumulative Distribution Analysis ===\n")
  cat(paste("Watershed:", watershed, "\n"))
  cat(paste("Years:", paste(years, collapse = ", "), "\n"))
  cat(paste("Data points every", interval_days, "days\n"))
  
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
  
  ################################################################################
  # PROCESS EACH YEAR
  ################################################################################
  
  for (year in years) {
    cat(paste("\n--- Processing year", year, "---\n"))
    
    #--------------------------------------------------------------------------#
    # 1. LOAD DATA
    #--------------------------------------------------------------------------#
    
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
    
    #--------------------------------------------------------------------------#
    # 2. SETUP ASSIGNMENT PARAMETERS
    #--------------------------------------------------------------------------#
    
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, params$min_error)
    priors <- setup_watershed_priors(edges, params$min_stream_order, watershed, natal_data)
    
    #--------------------------------------------------------------------------#
    # 3. DETERMINE DOY RANGE AND CREATE 3-DAY INTERVALS
    #--------------------------------------------------------------------------#
    
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
    cat(paste("Breakpoints:", paste(head(doy_breakpoints, 10), collapse = ", "), "...\n"))
    
    #--------------------------------------------------------------------------#
    # 4. CALCULATE CUMULATIVE DISTRIBUTIONS AT EACH 3-DAY INTERVAL
    #--------------------------------------------------------------------------#
    
    year_cumulative_data <- data.frame()
    
    for (i in 1:length(doy_breakpoints)) {
      current_doy <- doy_breakpoints[i]
      
      # Filter data up to current DOY (cumulative)
      cumulative_data_subset <- natal_data %>%
        filter(DOY <= current_doy)
      
      if (nrow(cumulative_data_subset) == 0) {
        # If no data yet, create zero entries for all management units
        cat(paste("No data yet at DOY", current_doy, "- skipping\n"))
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
            interval_number = i,  # Track which 3-day interval this is
            cumulative_production = total_production,
            fish_count_cumulative = nrow(cumulative_data_subset),
            total_production_all_units = total_production_all
          ) %>%
          select(year, doy, interval_number, mgmt_river, cumulative_production, 
                 fish_count_cumulative, total_production_all_units, edge_count)
        
        year_cumulative_data <- bind_rows(year_cumulative_data, doy_data)
      }
    }
    
    #--------------------------------------------------------------------------#
    # 5. CALCULATE FINAL PRODUCTION AND CUMULATIVE PERCENTAGES FOR EACH MGMT UNIT
    #--------------------------------------------------------------------------#
    
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
      #--------------------------------------------------------------------------#
      # 6. CALCULATE SUMMARY STATISTICS FOR THIS YEAR
      #--------------------------------------------------------------------------#
      
      # For each management unit, find key timing statistics
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
      
      all_summary_stats[[as.character(year)]] <- year_summary
      all_cumulative_data[[as.character(year)]] <- year_cumulative_data
    }
    
    cat(paste("Completed year", year, "with", nrow(year_cumulative_data), "data points\n"))
  }
  
  ################################################################################
  # 6. COMBINE ALL YEARS AND CREATE PLOTS
  ################################################################################
  
  if (length(all_cumulative_data) == 0) {
    stop("No valid data processed for any year")
  }
  
  # Combine all cumulative data
  combined_cumulative_data <- bind_rows(all_cumulative_data)
  combined_summary_stats <- bind_rows(all_summary_stats)
  
  cat(paste("Combined data:", nrow(combined_cumulative_data), "total data points\n"))
  cat(paste("Management units found:", length(unique(combined_cumulative_data$mgmt_river)), "\n"))
  
  ################################################################################
  # 7. CREATE VISUALIZATION PLOTS
  ################################################################################
  
  plots <- create_cumulative_distribution_plots(combined_cumulative_data, combined_summary_stats, 
                                                watershed, interval_days)
  
  ################################################################################
  # 8. EXPORT DATA IF REQUESTED
  ################################################################################
  
  if (export_csv) {
    export_cumulative_data(combined_cumulative_data, combined_summary_stats, 
                           watershed, interval_days)
  }
  
  ################################################################################
  # 9. RETURN RESULTS
  ################################################################################
  
  return(list(
    cumulative_data = combined_cumulative_data,
    summary_stats = combined_summary_stats,
    plots = plots,
    parameters = list(
      watershed = watershed,
      years = years,
      interval_days = interval_days,
      params = params
    )
  ))
}

#' Create comprehensive cumulative distribution plots
create_cumulative_distribution_plots <- function(cumulative_data, summary_stats, 
                                                 watershed, interval_days) {
  
  cat("Creating visualization plots...\n")
  
  # Create output directory
  plot_dir <- here("Analysis_Results/Cumulative_Distribution/Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # WATERSHED GRADIENT COLORING - RED → ORANGE → BLUE
  # Order from highest to lowest in watershed
  desired_order <- c(
    "N. Fork Kusko",
    "E. Fork Kuskokwim", 
    "S. Fork Kusko", 
    "Takotna and Nixon Fork",
    "Big River",
    "Upper Kusko Main",
    "Tatlawiksuk",
    "Kwethluk",
    "Stony",
    "Swift",
    "Holitna and Hoholitna", 
    "George",
    "Oskakawlik",
    "Middle Kusko Main",
    "Holokuk",
    "Aniak",
    "Tuluksak",
    "Kisaralik",
    "Hoholitna",
    "Johnson",
    "Lower Kusko"
  )
  
  # Create red → orange → blue gradient function
  create_watershed_colors <- function(watershed_order, mgmt_rivers_present) {
    # Filter to only management units present in data
    filtered_order <- watershed_order[watershed_order %in% mgmt_rivers_present]
    n_watersheds <- length(filtered_order)
    
    # Create gradient from red → orange → blue
    color_gradient <- colorRampPalette(c(
      "#8B0000",  # Dark red (highest watersheds)
      "#CC0000",  # Red
      "#FF0000",  # Bright red
      "#FF4500",  # Orange red
      "#FF8C00",  # Dark orange
      "#FFA500",  # Orange
      "#FFB347",  # Light orange
      "#87CEEB",  # Sky blue
      "#4682B4",  # Steel blue
      "#1E90FF",  # Dodger blue
      "#0000FF",  # Blue
      "#000080"   # Navy blue (lowest watersheds)
    ))(n_watersheds)
    
    # Create named vector
    names(color_gradient) <- filtered_order
    return(color_gradient)
  }
  
  # Generate colors for management units present in data
  mgmt_rivers <- unique(cumulative_data$mgmt_river)
  watershed_colors <- create_watershed_colors(desired_order, mgmt_rivers)
  
  # Create ordered factor for proper legend ordering (red to blue)
  mgmt_rivers_ordered <- desired_order[desired_order %in% mgmt_rivers]
  cumulative_data$mgmt_river <- factor(cumulative_data$mgmt_river, levels = mgmt_rivers_ordered)
  
  cat(paste("Applied red → orange → blue gradient to", length(watershed_colors), "management units\n"))
  cat("Legend will be ordered from red (highest) to blue (lowest)\n")
  
  #--------------------------------------------------------------------------#
  # PLOT 1: INDIVIDUAL PLOTS FOR EACH YEAR - ENHANCED FOR TIMING ANALYSIS
  #--------------------------------------------------------------------------#
  
  individual_plots <- list()
  years <- sort(unique(cumulative_data$year))
  
  cat(paste("Creating individual plots for", length(years), "years...\n"))
  
  for (year in years) {
    year_data <- cumulative_data %>% filter(year == !!year)
    
    # Load natal data for CPUE histogram
    natal_data <- load_natal_data(year, watershed)
    
    # Create CPUE histogram for this year
    cpue_histogram <- create_cpue_histogram(natal_data, year, watershed)
    
    # Create enhanced plot with timing grid and reference lines
    p_year <- ggplot(year_data, aes(x = doy, y = cumulative_percent_of_unit, color = mgmt_river)) +
      # Add reference lines for key timing milestones
      geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.7) +
      geom_hline(yintercept = c(10, 90), color = "gray80", linetype = "dotted", alpha = 0.5) +
      
      # Main data
      geom_line(linewidth = 2, alpha = 0.9) +
      geom_point(size = 3, alpha = 0.8) +
      
      # Color and scales
      scale_color_manual(values = watershed_colors, 
                         breaks = mgmt_rivers_ordered,  # Order legend red to blue
                         guide = guide_legend(override.aes = list(linewidth = 3))) +
      scale_y_continuous(
        labels = function(x) paste0(round(x, 1), "%"),
        limits = c(0, 100),
        breaks = c(0, 10, 25, 50, 75, 90, 100),
        minor_breaks = seq(0, 100, 5)
      ) +
      scale_x_continuous(
        limits = c(min(year_data$doy), max(year_data$doy)),
        breaks = scales::pretty_breaks(n = 8),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      
      # Labels and theme
      labs(
        title = paste("Salmon Run Timing Progress -", year),
        subtitle = paste(watershed, "Watershed | DOY", min(year_data$doy), "to", max(year_data$doy), "| Red (headwaters) → Orange → Blue (mouth)"),
        x = paste("Day of Year (DOY) - Range:", min(year_data$doy), "to", max(year_data$doy)),
        y = "Cumulative Percent of Management Unit's Total Production",
        color = "Management Unit\n(Watershed Position)",
        caption = paste("Data points every", interval_days, "days | Dashed lines: 25%, 50%, 75% milestones | Color gradient: Red=upstream, Blue=downstream")
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "gray20"),
        plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 10, hjust = 0.5, color = "gray50"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 11)
      )
    
    # Save ONLY the combined plot (this is what you want)
    combined_filename <- paste0("run_timing_progress_", year, "_", watershed, ".png")
    
    # Create and save combined plot
    png(file.path(plot_dir, combined_filename), width = 14, height = 11, units = "in", res = 300, bg = "white")
    grid.arrange(p_year, cpue_histogram, nrow = 2, heights = c(3, 1))
    dev.off()
    
    individual_plots[[as.character(year)]] <- p_year
    cat(paste("✓ Created combined timing plot for", year, "| DOY range:", min(year_data$doy), "to", max(year_data$doy), "| CPUE axis: June 1 - July 30 | File:", combined_filename, "\n"))
  }
  
  # Order data for consistent legend ordering
  cumulative_data$mgmt_river <- factor(cumulative_data$mgmt_river, levels = mgmt_rivers_ordered)
  
  # Set consistent DOY range for overview plot: May 28 (DOY 148) to July 29 (DOY 210)
  may_28_doy <- 148  # May 28
  july_29_doy <- 210  # July 29
  
  p1_faceted <- ggplot(cumulative_data, aes(x = doy, y = cumulative_percent_of_unit, color = mgmt_river)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.7) +
    facet_wrap(~year, ncol = 1) +  # Stack vertically (1 column)
    scale_color_manual(values = watershed_colors,
                       breaks = mgmt_rivers_ordered,  # Order legend red to blue
                       guide = guide_legend(override.aes = list(linewidth = 2))) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, 100),
      breaks = seq(0, 100, 25)
    ) +
    scale_x_continuous(
      limits = c(may_28_doy, july_29_doy),  # Fixed limits for all years
      breaks = seq(may_28_doy, july_29_doy, by = 14),  # Every 2 weeks
      labels = function(x) {
        # Convert DOY to date
        date_obj <- as.Date(x - 1, origin = "2024-01-01")
        format(date_obj, "%b %d")
      },
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = paste("Run Timing Progress Overview -", watershed, "Watershed"),
      subtitle = paste("May 28 to July 29 | Red (headwaters) → Orange → Blue (mouth)"),
      x = "Day of Year (Date) - Fixed Range: May 28 to July 29",
      y = "Cumulative Percent of Unit's Total Production",
      color = "Management Unit\n(Watershed Position)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )
  
  ggsave(file.path(plot_dir, paste0("cumulative_progress_overview_", watershed, ".png")), 
         p1_faceted, width = 12, height = 14, dpi = 300, bg = "white")  # Increased height for stacked layout
  
  #--------------------------------------------------------------------------#
  # PLOT 2: AVERAGE CUMULATIVE CURVES - ALSO WITH DATA-SCALED X-AXIS
  #--------------------------------------------------------------------------#
  
  avg_cumulative <- cumulative_data %>%
    group_by(mgmt_river, doy) %>%
    reframe(
      avg_cumulative_percent = mean(cumulative_percent_of_unit, na.rm = TRUE),
      sd_cumulative_percent = sd(cumulative_percent_of_unit, na.rm = TRUE),
      n_years = n()
    ) %>%
    mutate(
      se_cumulative_percent = sd_cumulative_percent / sqrt(n_years),
      ci_lower = pmax(0, avg_cumulative_percent - 1.96 * se_cumulative_percent),
      ci_upper = pmin(100, avg_cumulative_percent + 1.96 * se_cumulative_percent)
    )
  
  # Calculate DOY range for average plot
  avg_doy_range <- range(avg_cumulative$doy)
  
  p2 <- ggplot(avg_cumulative, aes(x = doy, y = avg_cumulative_percent, color = mgmt_river)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = mgmt_river), 
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.5, alpha = 0.9) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = watershed_colors,
                       breaks = mgmt_rivers_ordered,  # Order legend red to blue
                       guide = guide_legend(override.aes = list(linewidth = 2))) +
    scale_fill_manual(values = watershed_colors, guide = "none") +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, 100),
      breaks = seq(0, 100, 20)
    ) +
    scale_x_continuous(
      limits = avg_doy_range,
      breaks = scales::pretty_breaks(n = 8),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = paste("Average Run Timing Progress -", watershed, "Watershed"),
      subtitle = paste("Multi-year averages | DOY", avg_doy_range[1], "to", avg_doy_range[2], "| Red→Orange→Blue: upstream to downstream"),
      x = paste("Day of Year (DOY) - Data Range:", avg_doy_range[1], "to", avg_doy_range[2]),
      y = "Average Cumulative Percent of Unit's Total Production",
      color = "Management Unit\n(Watershed Position)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5)
    )
  
  ggsave(file.path(plot_dir, paste0("average_cumulative_progress_", watershed, ".png")), 
         p2, width = 12, height = 8, dpi = 300, bg = "white")
  
  #--------------------------------------------------------------------------#
  # PLOT 3: TIMING SUMMARY STATISTICS
  #--------------------------------------------------------------------------#
  
  # Calculate average timing stats across years
  avg_timing <- summary_stats %>%
    group_by(mgmt_river) %>%
    reframe(
      avg_median_timing = mean(median_timing, na.rm = TRUE),
      avg_timing_spread = mean(timing_spread, na.rm = TRUE),
      avg_final_percent = mean(final_overall_percent, na.rm = TRUE)
    ) %>%
    arrange(avg_median_timing)
  
  p3 <- ggplot(avg_timing, aes(x = reorder(mgmt_river, avg_median_timing))) +
    geom_col(aes(y = avg_median_timing, fill = mgmt_river), alpha = 0.8) +
    scale_fill_manual(values = watershed_colors, guide = "none") +
    coord_flip() +
    labs(
      title = paste("Average Median Timing by Management Unit -", watershed),
      subtitle = "Earlier timing = earlier in the run | Red (headwaters) → Blue (mouth)",
      x = "Management Unit (Watershed Position)",
      y = "Average Median Timing (DOY)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(plot_dir, paste0("timing_summary_", watershed, ".png")), 
         p3, width = 10, height = 8, dpi = 300, bg = "white")
  
  #--------------------------------------------------------------------------#
  # PLOT 4: TIMING SPREAD VS CONTRIBUTION
  #--------------------------------------------------------------------------#
  
  p4 <- ggplot(avg_timing, aes(x = avg_final_percent, y = avg_timing_spread, 
                               color = mgmt_river, size = avg_final_percent)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = watershed_colors,
                       breaks = mgmt_rivers_ordered,  # Order legend red to blue
                       guide = guide_legend(override.aes = list(linewidth = 2))) +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    scale_x_continuous(labels = function(x) paste0(round(x, 1), "%")) +
    labs(
      title = paste("Timing Spread vs Total Contribution -", watershed),
      subtitle = "Size = contribution size | Red (headwaters) → Orange → Blue (mouth)",
      x = "Average Total Contribution (%)",
      y = "Average Timing Spread (Days: 90th - 10th percentile)",
      color = "Management Unit\n(Watershed Position)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(plot_dir, paste0("timing_spread_vs_contribution_", watershed, ".png")), 
         p4, width = 12, height = 8, dpi = 300, bg = "white")
  
  cat(paste("Saved", length(individual_plots), "individual year plots with RED→ORANGE→BLUE gradient +", "4 summary plots to:", plot_dir, "\n"))
  
  return(list(
    individual_year_plots = individual_plots,
    overview_faceted = p1_faceted,
    average_cumulative = p2,
    timing_summary = p3,
    spread_vs_contribution = p4
  ))
}

#' Export cumulative distribution data to CSV files
export_cumulative_data <- function(cumulative_data, summary_stats, watershed, interval_days) {
  
  cat("Exporting data to CSV files...\n")
  
  # Create output directory
  csv_dir <- here("Analysis_Results/Cumulative_Distribution/CSV")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  #--------------------------------------------------------------------------#
  # EXPORT 1: FULL CUMULATIVE DATA (LONG FORMAT)
  #--------------------------------------------------------------------------#
  
  cumulative_clean <- cumulative_data %>%
    select(year, doy, interval_number, mgmt_river, cumulative_production, 
           cumulative_percent_of_unit, overall_percent, fish_count_cumulative) %>%
    arrange(year, mgmt_river, doy)
  
  cumulative_filepath <- file.path(csv_dir, paste0("cumulative_distribution_data_", 
                                                   watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(cumulative_clean, cumulative_filepath, row.names = FALSE)
  cat(paste("Exported full cumulative data to:", basename(cumulative_filepath), "\n"))
  
  #--------------------------------------------------------------------------#
  # EXPORT 2: SUMMARY STATISTICS
  #--------------------------------------------------------------------------#
  
  summary_filepath <- file.path(csv_dir, paste0("timing_summary_stats_", 
                                                watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(summary_stats, summary_filepath, row.names = FALSE)
  cat(paste("Exported summary statistics to:", basename(summary_filepath), "\n"))
  
  #--------------------------------------------------------------------------#
  # EXPORT 3: WIDE FORMAT FOR EASY PLOTTING
  #--------------------------------------------------------------------------#
  
  cumulative_wide <- cumulative_data %>%
    select(year, doy, mgmt_river, cumulative_percent_of_unit) %>%
    pivot_wider(names_from = mgmt_river, values_from = cumulative_percent_of_unit, 
                values_fill = 0, names_prefix = "mgmt_") %>%
    arrange(year, doy)
  
  wide_filepath <- file.path(csv_dir, paste0("cumulative_distribution_wide_", 
                                             watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(cumulative_wide, wide_filepath, row.names = FALSE)
  cat(paste("Exported wide format data to:", basename(wide_filepath), "\n"))
  
  #--------------------------------------------------------------------------#
  # EXPORT 4: AVERAGE ACROSS YEARS
  #--------------------------------------------------------------------------#
  
  avg_cumulative_export <- cumulative_data %>%
    group_by(mgmt_river, doy) %>%
    reframe(
      avg_cumulative_percent = mean(cumulative_percent_of_unit, na.rm = TRUE),
      sd_cumulative_percent = sd(cumulative_percent_of_unit, na.rm = TRUE),
      n_years = n(),
      min_cumulative_percent = min(cumulative_percent_of_unit, na.rm = TRUE),
      max_cumulative_percent = max(cumulative_percent_of_unit, na.rm = TRUE)
    ) %>%
    arrange(mgmt_river, doy)
  
  avg_filepath <- file.path(csv_dir, paste0("average_cumulative_distribution_", 
                                            watershed, "_", interval_days, "day_intervals.csv"))
  write.csv(avg_cumulative_export, avg_filepath, row.names = FALSE)
  cat(paste("Exported average cumulative data to:", basename(avg_filepath), "\n"))
  
  cat(paste("All CSV exports completed in:", csv_dir, "\n"))
}

################################################################################
# MAIN EXECUTION - RUNS AUTOMATICALLY WHEN SCRIPT IS EXECUTED
################################################################################

cat("################################################################################\n")
cat("# EXECUTING CUMULATIVE DISTRIBUTION ANALYSIS\n")
cat("################################################################################\n")

# Run the analysis
cat("Starting Cumulative Distribution Analysis...\n")
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
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat(paste("Total data points:", nrow(results$cumulative_data), "\n"))
  cat(paste("Management units analyzed:", length(unique(results$cumulative_data$mgmt_river)), "\n"))
  cat(paste("DOY range:", min(results$cumulative_data$doy), "to", max(results$cumulative_data$doy), "\n"))
  cat(paste("Data points created every", results$parameters$interval_days, "days\n"))
  
  # Verify 3-day intervals
  doy_intervals <- results$cumulative_data %>%
    filter(mgmt_river == first(mgmt_river)) %>%  # Pick one management unit
    filter(year == first(year)) %>%              # Pick one year
    arrange(doy) %>%
    mutate(interval_gap = doy - lag(doy, default = first(doy))) %>%
    filter(row_number() > 1)  # Remove first row (no lag)
  
  cat("\n=== INTERVAL VERIFICATION ===\n")
  cat("Gaps between consecutive data points (should be 3 days):\n")
  print(table(doy_intervals$interval_gap))
  
  # Display sample results
  cat("\n=== SAMPLE RESULTS ===\n")
  cat("First few rows of cumulative data:\n")
  print(head(results$cumulative_data, 10))
  
  cat("\n=== SUMMARY STATISTICS SAMPLE ===\n")
  cat("Sample timing statistics:\n")
  print(head(results$summary_stats, 5))
  
  cat("\n=== FILES CREATED ===\n")
  cat("📊 INDIVIDUAL YEAR TIMING PLOTS:\n")
  years_processed <- sort(unique(results$cumulative_data$year))
  for (year in years_processed) {
    filename <- paste0("run_timing_progress_", year, "_", results$parameters$watershed, ".png")
    cat(paste("   ", year, ":", filename, "\n"))
  }
  cat(paste("\n📁 Plot Location:", here("Analysis_Results/Cumulative_Distribution/Plots/"), "\n"))
  cat("📁 CSV files saved to: Analysis_Results/Cumulative_Distribution/CSV/\n")
  cat("\n🎯 USE THESE PLOTS TO SEE WHEN EACH POPULATION ARRIVES!\n")
  cat("   🎨 COLOR GRADIENT: Red (headwaters) → Orange (middle) → Blue (mouth)\n")
  cat("   - Each curve shows 0% → 100% progression for that management unit\n")
  cat("   - Steep curves = populations arrive quickly\n") 
  cat("   - Gradual curves = populations spread out over time\n")
  cat("   - Compare timing between management units within each year\n")
  cat("   - Compare same management unit across different years\n")
  cat("   - Red units = upstream/headwater populations\n")
  cat("   - Blue units = downstream/mouth populations\n")
  cat("\nAnalysis completed successfully!\n")
  
} else {
  cat("Analysis failed. Please check error messages above.\n")
}

cat("################################################################################\n")
cat("# SCRIPT EXECUTION COMPLETE\n")
cat("################################################################################\n")