# Analysis_utils.R
# Utility functions for DFA calculation and analysis

#' Process DFA data for multiple years and quartiles
#' 
#' @param watershed Character, watershed name
#' @param years Numeric vector, years to process
#' @param sensitivity_threshold Numeric, sensitivity threshold for assignment
#' @param min_error Numeric, minimum error value
#' @param min_stream_order Numeric, minimum stream order
#' @return List containing combined results and time series data
process_dfa_data <- function(watershed, years, sensitivity_threshold = 0.7, 
                             min_error = 0.0006, min_stream_order = 3) {
  
  all_results <- list()
  
  for (year in years) {
    cat("Processing year", year, "\n")
    
    # Load data
    spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
    edges <- spatial_data$edges
    basin <- spatial_data$basin
    Huc <- spatial_data$Huc
    natal_data <- load_natal_data(year, watershed)
    
    # Calculate total year CPUE
    total_year_cpue <- sum(natal_data$dailyCPUEprop, na.rm = TRUE)
    
    # Divide into quartiles
    quartile_data <- divide_doy_quartiles(natal_data)
    quartile_subsets <- quartile_data$subsets
    
    # Set up assignment parameters
    pid_iso <- edges$iso_pred
    pid_isose <- edges$isose_pred
    error <- calculate_error(pid_isose, min_error)
    priors <- setup_watershed_priors(edges, min_stream_order, watershed, natal_data)
    
    # Process each quartile
    for (q in 1:4) {
      current_subset <- quartile_subsets[[q]]
      
      if (nrow(current_subset) == 0) next
      
      # Calculate quartile CPUE proportion
      quartile_cpue <- sum(current_subset$dailyCPUEprop, na.rm = TRUE)
      quartile_cpue_proportion <- quartile_cpue / total_year_cpue
      
      # Perform assignment
      assignment_matrix <- perform_assignment(
        current_subset, edges, watershed, priors, pid_iso, error, sensitivity_threshold
      )
      
      # Calculate basin assignments
      basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
      
      # Process HUC data
      huc_result <- process_huc_data(edges, basin, Huc, basin_assign_sum, 8)
      
      # Calculate proportions
      total_huc_production <- sum(huc_result$total_production, na.rm = TRUE)
      huc_result$huc_proportion_within_quartile <- huc_result$total_production / total_huc_production
      huc_result$contribution_to_total_year <- huc_result$huc_proportion_within_quartile * quartile_cpue_proportion
      
      # Store results
      result_data <- st_drop_geometry(huc_result) %>%
        select(Name, huc_proportion_within_quartile, contribution_to_total_year) %>%
        mutate(
          year = year,
          quartile = q,
          quartile_cpue_proportion = quartile_cpue_proportion
        )
      
      key <- paste(year, "Q", q, sep = "_")
      all_results[[key]] <- result_data
    }
  }
  
  # Combine all results
  combined_data <- bind_rows(all_results)
  
  # Create time series data (HUCs only)
  time_series <- combined_data %>%
    select(Name, year, quartile, contribution_to_total_year) %>%
    mutate(time_step = paste0(year, "_Q", quartile)) %>%
    select(Name, time_step, contribution_to_total_year)
  
  # Create total timeseries for plotting
  overall_run_timeseries <- combined_data %>%
    select(year, quartile, quartile_cpue_proportion) %>%
    distinct() %>%
    mutate(
      Name = "Total",
      time_step = paste0(year, "_Q", quartile),
      contribution_to_total_year = quartile_cpue_proportion
    ) %>%
    select(Name, time_step, contribution_to_total_year)
  
  time_series_with_total <- bind_rows(time_series, overall_run_timeseries)
  
  return(list(
    combined_data = combined_data,
    time_series = time_series,
    time_series_with_total = time_series_with_total
  ))
}

#' Create and save formatted data files
#' 
#' @param results List from process_dfa_data()
#' @param output_dir Character, output directory path
#' @return Data frame with z-normalized time series
save_dfa_files <- function(results, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save long format
  long_filepath <- file.path(output_dir, "kusko_dfa_long.csv")
  write.csv(results$combined_data, long_filepath, row.names = FALSE)
  
  # Create and save wide format
  time_series_wide <- results$time_series %>%
    pivot_wider(names_from = time_step, values_from = contribution_to_total_year, values_fill = 0)
  
  wide_filepath <- file.path(output_dir, "kusko_dfa_wide.csv")
  write.csv(time_series_wide, wide_filepath, row.names = FALSE)
  
  # Create z-normalized version
  huc_names <- time_series_wide$Name
  ts_values <- time_series_wide[, -1]
  z_normalized_ts <- t(scale(t(ts_values)))
  
  time_series_z_normalized <- data.frame(
    Name = huc_names,
    z_normalized_ts,
    stringsAsFactors = FALSE
  )
  
  # Remove the 5th row (as in original)
  time_series_z_normalized <- time_series_z_normalized[-5, ]
  
  z_normalized_filepath <- file.path(output_dir, "kusko_dfa_wide_znormalized.csv")
  write.csv(time_series_z_normalized, z_normalized_filepath, row.names = FALSE)
  
  cat("Files saved to:", output_dir, "\n")
  cat("Z-normalized timeseries saved to:", z_normalized_filepath, "\n")
  
  return(time_series_z_normalized)
}

#' Create DFA visualization plots
#' 
#' @param time_series_with_total Data frame with time series including total
#' @param time_series_z_normalized Data frame with z-normalized data
#' @return List of ggplot objects
create_dfa_plots <- function(time_series_with_total, time_series_z_normalized) {
  
  # Prepare plot data
  plot_data <- time_series_with_total %>%
    separate(time_step, into = c("year", "quartile"), sep = "_Q", remove = FALSE) %>%
    mutate(
      year_numeric = as.numeric(year),
      quartile_numeric = as.numeric(quartile),
      time_continuous = year_numeric + (quartile_numeric - 1) * 0.25
    )
  
  # Plot 1: Overall run timing
  plot_data_total <- plot_data %>% filter(Name == "Total")
  
  p1 <- ggplot(plot_data_total, aes(x = time_continuous, y = contribution_to_total_year * 100)) +
    geom_line(linewidth = 2, color = "black") +
    geom_point(size = 3, color = "black") +
    scale_x_continuous(breaks = seq(2017, 2021, 1)) +
    labs(
      title = "Overall Salmon Run Timing Pattern",
      x = "Year",
      y = "Proportion of Annual Run (%)"
    ) +
    theme_minimal()
  
  # Plot 2: Individual HUC contributions
  plot_data_hucs <- plot_data %>% filter(Name != "Total")
  
  p2 <- ggplot(plot_data_hucs, aes(x = time_continuous, y = contribution_to_total_year * 100, color = Name)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(plot_data_hucs$Name)))) +
    scale_x_continuous(breaks = seq(2017, 2021, 1)) +
    labs(
      title = "Individual HUC Contributions",
      x = "Year",
      y = "Contribution to Annual Run (%)",
      color = "HUC"
    ) +
    theme_minimal()
  
  # Plot 3: Z-normalized faceted plot
  z_normalized_long <- time_series_z_normalized %>%
    pivot_longer(cols = -Name, names_to = "time_step", values_to = "z_normalized_value") %>%
    separate(time_step, into = c("year", "quartile"), sep = "_Q", remove = FALSE) %>%
    mutate(
      year_numeric = as.numeric(year),
      quartile_numeric = as.numeric(quartile),
      time_continuous = year_numeric + (quartile_numeric - 1) * 0.25
    )
  
  p3 <- ggplot(z_normalized_long, aes(x = time_continuous, y = z_normalized_value)) +
    geom_line(linewidth = 1, color = "blue") +
    geom_point(size = 1.5, color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
    facet_wrap(~Name, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = seq(2017, 2021, 2)) +
    labs(
      title = "Z-Normalized HUC Timeseries",
      subtitle = "Each HUC standardized to mean=0, sd=1",
      x = "Year",
      y = "Z-normalized Contribution"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    )
  
  return(list(
    overall_timing = p1,
    huc_contributions = p2,
    z_normalized_facets = p3
  ))
}

#' Fit DFA model to z-normalized data
#' 
#' @param time_series_z_normalized Data frame with z-normalized time series
#' @param state_numbers Numeric, number of states for DFA model
#' @param model_params List, model parameters for MARSS
#' @param control_params List, control parameters for MARSS
#' @return MARSS model object
fit_dfa_model <- function(time_series_z_normalized, state_numbers = 5, 
                          model_params = list(m = 3, R = "diagonal and unequal"),
                          control_params = list(maxit = 50000)) {
  
  # Convert to matrix format
  dfa_data_matrix <- as.matrix(time_series_z_normalized[, -1])
  rownames(dfa_data_matrix) <- time_series_z_normalized$Name
  
  # Fit the DFA model
  model <- suppressMessages(suppressWarnings(
    MARSS(dfa_data_matrix,
          model = model_params,
          z.score = FALSE,  # Data already z-normalized
          form = "dfa", 
          control = control_params, 
          method = "BFGS")
  ))
  
  return(model)
}