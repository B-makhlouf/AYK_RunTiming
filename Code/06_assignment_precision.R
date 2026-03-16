################################################################################
# 06_ASSIGNMENT_PRECISION.R - INDIVIDUAL ASSIGNMENT PRECISION ANALYSIS
################################################################################
# PURPOSE: Calculate and visualize precision of individual fish assignments
# Shows what percentage of each fish's probability falls in its top management unit
# Integrates seamlessly with existing workflow
# Run 00_setup.R first
################################################################################

cat("=== ASSIGNMENT PRECISION ANALYSIS ===\n")

# Check if setup is loaded
if (!exists("CONFIG")) {
  stop("Please run 00_setup.R first to load parameters and functions")
}

library(dplyr)
library(ggplot2)
library(glue)

################################################################################
# CORE PRECISION CALCULATION FUNCTIONS
################################################################################

#' Calculate individual assignment precision within management units
#' 
#' @param assignment_matrix Matrix from perform_assignment (edges × fish)
#' @param edges SF object with stream edges including mgmt_river column
#' @return Data frame with precision metrics for each individual
calculate_individual_precision <- function(assignment_matrix, edges) {
  
  # Check for management river data
  if (!"mgmt_river" %in% colnames(edges)) {
    warning("mgmt_river column not found in edges data")
    return(NULL)
  }
  
  # Get management unit for each edge
  mgmt_units <- edges$mgmt_river
  
  # Filter to only managed edges
  managed_idx <- which(!is.na(mgmt_units) & mgmt_units != "")
  
  if (length(managed_idx) == 0) {
    warning("No managed edges found")
    return(NULL)
  }
  
  # Initialize results
  n_fish <- ncol(assignment_matrix)
  precision_results <- data.frame(
    fish_id = 1:n_fish,
    top_mgmt_unit = character(n_fish),
    pct_in_top_unit = numeric(n_fish),
    total_assignment = numeric(n_fish),
    n_units_assigned = integer(n_fish),
    stringsAsFactors = FALSE
  )
  
  # Process each fish
  for (i in 1:n_fish) {
    # Get assignments for this fish (only managed edges)
    fish_assignments <- assignment_matrix[managed_idx, i]
    fish_mgmt_units <- mgmt_units[managed_idx]
    
    # Remove NA and zero assignments
    valid_idx <- !is.na(fish_assignments) & fish_assignments > 0
    fish_assignments <- fish_assignments[valid_idx]
    fish_mgmt_units <- fish_mgmt_units[valid_idx]
    
    if (length(fish_assignments) == 0) {
      precision_results$top_mgmt_unit[i] <- NA
      precision_results$pct_in_top_unit[i] <- 0
      precision_results$total_assignment[i] <- 0
      precision_results$n_units_assigned[i] <- 0
      next
    }
    
    # Sum assignments by management unit
    mgmt_totals <- tapply(fish_assignments, fish_mgmt_units, sum, na.rm = TRUE)
    
    # Find top management unit
    top_unit <- names(mgmt_totals)[which.max(mgmt_totals)]
    top_unit_total <- max(mgmt_totals)
    
    # Calculate percentage in top unit
    total_assignment <- sum(fish_assignments, na.rm = TRUE)
    pct_in_top <- (top_unit_total / total_assignment) * 100
    
    # Store results
    precision_results$top_mgmt_unit[i] <- top_unit
    precision_results$pct_in_top_unit[i] <- pct_in_top
    precision_results$total_assignment[i] <- total_assignment
    precision_results$n_units_assigned[i] <- length(mgmt_totals)
  }
  
  return(precision_results)
}

################################################################################
# VISUALIZATION FUNCTIONS
################################################################################

#' Create histogram of individual assignment precision
#' 
#' @param precision_data Data frame from calculate_individual_precision()
#' @param year Year being analyzed
#' @param watershed Watershed name
#' @param quartile Optional quartile label (e.g., "Q1")
#' @param binwidth Width of histogram bins (default 5%)
#' @return ggplot object
plot_precision_histogram <- function(precision_data, year, watershed, 
                                     quartile = NULL, binwidth = 5) {
  
  # Remove fish with no assignments
  plot_data <- precision_data %>%
    filter(!is.na(pct_in_top_unit) & total_assignment > 0)
  
  if (nrow(plot_data) == 0) {
    warning("No valid data to plot")
    return(NULL)
  }
  
  # Calculate summary statistics
  mean_pct <- mean(plot_data$pct_in_top_unit, na.rm = TRUE)
  median_pct <- median(plot_data$pct_in_top_unit, na.rm = TRUE)
  n_fish <- nrow(plot_data)
  
  # Build title
  title_parts <- c(watershed, year)
  if (!is.null(quartile)) title_parts <- c(title_parts, quartile)
  
  main_title <- paste("Individual Assignment Precision:", paste(title_parts, collapse = " "))
  subtitle <- sprintf("n = %d fish | Mean = %.1f%% | Median = %.1f%%", 
                      n_fish, mean_pct, median_pct)
  
  # Create histogram
  p <- ggplot(plot_data, aes(x = pct_in_top_unit)) +
    geom_histogram(binwidth = binwidth, fill = "steelblue", 
                   color = "white", alpha = 0.8) +
    geom_vline(aes(xintercept = mean_pct), 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = median_pct), 
               color = "darkgreen", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    labs(
      title = main_title,
      subtitle = subtitle,
      x = "Percentage of Assignment in Top Management Unit (%)",
      y = "Number of Fish",
      caption = "Red line = mean | Green line = median\nAfter applying priors and 70% threshold"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create combined faceted histogram for all years
#' 
#' @param all_precision_data Combined data frame with year column
#' @param watershed Watershed name
#' @param binwidth Width of histogram bins
#' @return ggplot object
plot_precision_by_year <- function(all_precision_data, watershed, binwidth = 5) {
  
  # Remove invalid data
  plot_data <- all_precision_data %>%
    filter(!is.na(pct_in_top_unit) & total_assignment > 0)
  
  if (nrow(plot_data) == 0) {
    warning("No valid data to plot")
    return(NULL)
  }
  
  # Calculate summary stats by year
  summary_stats <- plot_data %>%
    group_by(year) %>%
    summarise(
      mean_pct = mean(pct_in_top_unit, na.rm = TRUE),
      median_pct = median(pct_in_top_unit, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # Create faceted histogram
  p <- ggplot(plot_data, aes(x = pct_in_top_unit)) +
    geom_histogram(binwidth = binwidth, fill = "steelblue", 
                   color = "white", alpha = 0.8) +
    geom_vline(data = summary_stats, aes(xintercept = mean_pct), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_vline(data = summary_stats, aes(xintercept = median_pct), 
               color = "darkgreen", linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ year, scales = "free_y", ncol = 2) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(
      title = paste("Individual Assignment Precision by Year -", watershed),
      x = "Percentage of Assignment in Top Management Unit (%)",
      y = "Number of Fish",
      caption = "Red line = mean | Green line = median"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

################################################################################
# SUMMARY STATISTICS FUNCTIONS
################################################################################

#' Generate summary table of precision metrics
#' 
#' @param precision_data Data frame from calculate_individual_precision()
#' @return Data frame with summary statistics
summarize_precision <- function(precision_data) {
  
  # Filter valid data
  valid_data <- precision_data %>%
    filter(!is.na(pct_in_top_unit) & total_assignment > 0)
  
  if (nrow(valid_data) == 0) {
    return(NULL)
  }
  
  # Calculate summary statistics
  summary <- valid_data %>%
    summarise(
      n_fish = n(),
      mean_pct = mean(pct_in_top_unit),
      median_pct = median(pct_in_top_unit),
      sd_pct = sd(pct_in_top_unit),
      min_pct = min(pct_in_top_unit),
      max_pct = max(pct_in_top_unit),
      q25_pct = quantile(pct_in_top_unit, 0.25),
      q75_pct = quantile(pct_in_top_unit, 0.75),
      pct_above_90 = sum(pct_in_top_unit > 90) / n() * 100,
      pct_above_80 = sum(pct_in_top_unit > 80) / n() * 100,
      pct_above_70 = sum(pct_in_top_unit > 70) / n() * 100,
      pct_below_50 = sum(pct_in_top_unit < 50) / n() * 100,
      mean_n_units = mean(n_units_assigned)
    )
  
  return(summary)
}

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

#' Run complete assignment precision analysis
#' Integrates with existing workflow and processes all years
#' 
#' @param years Vector of years to analyze
#' @param watersheds Vector of watersheds to analyze
#' @param export_csv Whether to export CSV files (default TRUE)
#' @return List containing all precision data and summaries
run_precision_analysis <- function(years = CONFIG$years,
                                   watersheds = CONFIG$watersheds,
                                   export_csv = TRUE) {
  
  cat("\n=== ASSIGNMENT PRECISION ANALYSIS ===\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n\n")
  
  # Create output directory
  output_dir <- file.path(PATHS$output_dir, "Assignment_Precision")
  fig_dir <- file.path(output_dir, "Figures")
  csv_dir <- file.path(output_dir, "CSV")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_precision_data <- list()
  all_summaries <- list()
  
  for (watershed in watersheds) {
    params <- WATERSHED_PARAMS[[watershed]]
    watershed_precision <- data.frame()
    
    for (year in years) {
      cat(glue("Processing precision for {watershed} {year}...\n"))
      
      # Load data
      spatial_data <- load_spatial_data(watershed)
      natal_data <- load_natal_data(year, watershed)
      
      cat(glue("  Loaded {nrow(natal_data)} fish observations\n"))
      
      # Check for management units
      if (!"mgmt_river" %in% colnames(spatial_data$edges)) {
        warning("Management river data not found. Skipping year ", year)
        next
      }
      
      # Setup assignment parameters
      pid_iso <- spatial_data$edges$iso_pred
      pid_isose <- spatial_data$edges$isose_pred
      error <- calculate_error(pid_isose, params$min_error)
      priors <- setup_priors(spatial_data$edges, watershed, natal_data)
      
      # Perform assignment
      cat("  Performing assignment...\n")
      assignment_matrix <- perform_assignment(
        natal_data, spatial_data$edges, watershed, priors,
        pid_iso, error, params$sensitivity_threshold
      )
      
      # Calculate precision
      cat("  Calculating individual precision...\n")
      precision_results <- calculate_individual_precision(
        assignment_matrix, spatial_data$edges
      )
      
      if (is.null(precision_results)) {
        warning("  Failed to calculate precision for ", year)
        next
      }
      
      # Add metadata
      precision_results$year <- year
      precision_results$watershed <- watershed
      
      # Store results
      all_precision_data[[glue("{watershed}_{year}")]] <- precision_results
      watershed_precision <- bind_rows(watershed_precision, precision_results)
      
      # Create individual year histogram
      cat("  Creating histogram...\n")
      p <- plot_precision_histogram(precision_results, year, watershed)
      
      if (!is.null(p)) {
        fig_path <- file.path(fig_dir, 
                              glue("{watershed}_{year}_precision_histogram.png"))
        ggsave(fig_path, p, width = 10, height = 6, dpi = 300)
      }
      
      # Calculate summary statistics
      summary <- summarize_precision(precision_results)
      if (!is.null(summary)) {
        summary$year <- year
        summary$watershed <- watershed
        all_summaries[[glue("{watershed}_{year}")]] <- summary
      }
      
      cat(glue("  ✓ Completed {watershed} {year}\n\n"))
    }
    
    # Create combined plot for all years in this watershed
    if (nrow(watershed_precision) > 0) {
      cat(glue("Creating combined plot for {watershed}...\n"))
      p_combined <- plot_precision_by_year(watershed_precision, watershed)
      
      if (!is.null(p_combined)) {
        fig_path <- file.path(fig_dir, 
                              glue("{watershed}_all_years_precision.png"))
        ggsave(fig_path, p_combined, width = 12, height = 8, dpi = 300)
      }
    }
  }
  
  # Combine all results
  if (length(all_precision_data) > 0) {
    combined_precision <- bind_rows(all_precision_data)
    combined_summaries <- bind_rows(all_summaries)
    
    # Export CSVs if requested
    if (export_csv) {
      cat("\nExporting CSV files...\n")
      
      # Individual fish data
      precision_file <- file.path(csv_dir, "individual_fish_precision.csv")
      write.csv(combined_precision, precision_file, row.names = FALSE)
      cat(glue("  ✓ Individual data: {precision_file}\n"))
      
      # Summary statistics
      summary_file <- file.path(csv_dir, "precision_summary_statistics.csv")
      write.csv(combined_summaries, summary_file, row.names = FALSE)
      cat(glue("  ✓ Summary stats: {summary_file}\n"))
      
      # Create detailed summary by year
      year_summary <- combined_precision %>%
        filter(!is.na(pct_in_top_unit) & total_assignment > 0) %>%
        group_by(watershed, year) %>%
        summarise(
          n_fish = n(),
          mean_precision = mean(pct_in_top_unit),
          median_precision = median(pct_in_top_unit),
          sd_precision = sd(pct_in_top_unit),
          pct_high_precision = sum(pct_in_top_unit > 80) / n() * 100,
          pct_low_precision = sum(pct_in_top_unit < 50) / n() * 100,
          mean_units_per_fish = mean(n_units_assigned),
          .groups = "drop"
        )
      
      year_summary_file <- file.path(csv_dir, "precision_by_year.csv")
      write.csv(year_summary, year_summary_file, row.names = FALSE)
      cat(glue("  ✓ Year summary: {year_summary_file}\n"))
    }
    
    # Print overall summary
    cat("\n=== OVERALL PRECISION SUMMARY ===\n")
    overall_summary <- combined_precision %>%
      filter(!is.na(pct_in_top_unit) & total_assignment > 0) %>%
      summarise(
        total_fish = n(),
        mean_precision = mean(pct_in_top_unit),
        median_precision = median(pct_in_top_unit),
        pct_high_precision_80 = sum(pct_in_top_unit > 80) / n() * 100,
        pct_high_precision_90 = sum(pct_in_top_unit > 90) / n() * 100,
        pct_low_precision_50 = sum(pct_in_top_unit < 50) / n() * 100
      )
    
    cat(glue("Total fish analyzed: {overall_summary$total_fish}\n"))
    cat(glue("Mean precision: {round(overall_summary$mean_precision, 1)}%\n"))
    cat(glue("Median precision: {round(overall_summary$median_precision, 1)}%\n"))
    cat(glue("Fish with >80% precision: {round(overall_summary$pct_high_precision_80, 1)}%\n"))
    cat(glue("Fish with >90% precision: {round(overall_summary$pct_high_precision_90, 1)}%\n"))
    cat(glue("Fish with <50% precision: {round(overall_summary$pct_low_precision_50, 1)}%\n"))
    
    cat(glue("\n✓ Analysis complete. Results saved to:\n  {output_dir}\n\n"))
    
    return(list(
      precision_data = combined_precision,
      summaries = combined_summaries,
      overall = overall_summary
    ))
  } else {
    warning("No precision data was calculated")
    return(NULL)
  }
}

################################################################################
# EXECUTION SECTION
################################################################################

if (interactive() || !exists(".precision_script_executed")) {
  cat("Assignment precision analysis functions loaded.\n")
  cat("Available function:\n")
  cat("  run_precision_analysis()  # Calculate precision for all years\n\n")
  
  # Uncomment to run automatically
  # results <- run_precision_analysis()
  
  .precision_script_executed <- TRUE
}

cat("✓ Assignment precision analysis functions loaded\n")