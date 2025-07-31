# =============================================================================
# COMPREHENSIVE UTILITIES FOR AYK RUN TIMING ANALYSIS
# 
# All utility functions for spatial data, assignment, validation, and plotting
# Shared across all three analyses: quarterly timing, cumulative distribution, 
# and front-end closure
#
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
  library(grid)
  library(gridExtra)
})

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load spatial data for watershed analysis
#' 
#' @param watershed Character: "Kusko" or "Yukon"  
#' @param min_stream_order Numeric: minimum stream order to include
#' @param edges_path Character: full path to edges shapefile
#' @param basin_path Character: full path to basin shapefile
#' @return List with edges and basin spatial objects
load_spatial_data <- function(watershed, min_stream_order = 3, 
                              edges_path = NULL, basin_path = NULL) {
  
  if (is.null(edges_path) || is.null(basin_path)) {
    stop("Please provide full paths to edges_path and basin_path shapefiles")
  }
  
  if (!file.exists(edges_path)) {
    stop("Edges shapefile not found at: ", edges_path)
  }
  if (!file.exists(basin_path)) {
    stop("Basin shapefile not found at: ", basin_path)
  }
  
  cat("Loading spatial data for", watershed, "watershed...\n")
  edges <- st_read(edges_path, quiet = TRUE)
  basin <- st_read(basin_path, quiet = TRUE)
  
  # Transform and filter data
  edges <- st_transform(edges, st_crs(basin))
  edges <- edges[edges$Str_Order >= min_stream_order, ]
  
  cat("Loaded", nrow(edges), "stream edges with order >=", min_stream_order, "\n")
  
  return(list(edges = edges, basin = basin))
}

#' Load natal origins data for specific year and watershed
#' 
#' @param year Numeric: year to load (e.g., 2017)
#' @param watershed Character: "Kusko" or "Yukon"
#' @param data_dir Character: path to natal origins data directory
#' @return Data frame with natal origins data
load_natal_data <- function(year, watershed, data_dir = "data/natal_origins") {
  
  file_path <- file.path(data_dir, paste0(year, "_", watershed, "_Natal_Origins_Genetics_CPUE.csv"))
  
  if (!file.exists(file_path)) {
    stop("Natal data file not found: ", file_path)
  }
  
  natal_data <- read.csv(file_path)
  
  # Clean data by removing rows with NA in key columns
  if (watershed == "Yukon") {
    clean_data <- natal_data[!is.na(natal_data$Lower) & 
                               !is.na(natal_data$natal_iso) & 
                               !is.na(natal_data$dailyCPUEprop), ]
  } else {
    clean_data <- natal_data[!is.na(natal_data$natal_iso) & 
                               !is.na(natal_data$dailyCPUEprop), ]
  }
  
  cat("Loaded", nrow(clean_data), "clean natal origin records for", year, watershed, "\n")
  return(clean_data)
}

# =============================================================================
# MANAGEMENT UNIT VALIDATION FUNCTIONS
# =============================================================================

#' Get standardized watershed ordering (upstream to downstream)
#' 
#' @return Character vector of all 21 management unit names in proper order
get_watershed_order <- function() {
  c(
    "N. Fork Kusko",
    "E. Fork Kuskokwim River", 
    "S. Fork Kusko",
    "Upper Kusko Main",
    "Big River",
    "Takotna and Nixon Fork",
    "Tatlawiksuk",
    "Swift",
    "Stony", 
    "Holitna",
    "Hoholitna",
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
}

#' Validate that all 21 management units are present and properly named
#' 
#' @param data Data frame or vector containing management unit names
#' @param mgmt_col Character: name of management unit column (if data frame)
#' @param analysis_name Character: name of current analysis for error messages
#' @return List with validation results and corrected data
validate_management_units <- function(data, mgmt_col = "mgmt_river", analysis_name = "Analysis") {
  
  expected_units <- get_watershed_order()
  n_expected <- length(expected_units)
  
  # Extract management units from data
  if (is.data.frame(data)) {
    if (!mgmt_col %in% colnames(data)) {
      stop("Column '", mgmt_col, "' not found in data for ", analysis_name)
    }
    actual_units <- unique(data[[mgmt_col]][!is.na(data[[mgmt_col]]) & data[[mgmt_col]] != ""])
  } else {
    actual_units <- unique(data[!is.na(data) & data != ""])
  }
  
  # Name variations that should be corrected
  name_corrections <- list(
    "E. Fork Kuskokwim" = "E. Fork Kuskokwim River",
    "N Fork Kusko" = "N. Fork Kusko",
    "S Fork Kusko" = "S. Fork Kusko"
  )
  
  # Apply name corrections
  corrected_data <- data
  corrections_made <- character(0)
  
  for (wrong_name in names(name_corrections)) {
    correct_name <- name_corrections[[wrong_name]]
    if (is.data.frame(corrected_data)) {
      if (wrong_name %in% corrected_data[[mgmt_col]]) {
        corrected_data[[mgmt_col]][corrected_data[[mgmt_col]] == wrong_name] <- correct_name
        corrections_made <- c(corrections_made, paste(wrong_name, "→", correct_name))
      }
    } else {
      if (wrong_name %in% corrected_data) {
        corrected_data[corrected_data == wrong_name] <- correct_name
        corrections_made <- c(corrections_made, paste(wrong_name, "→", correct_name))
      }
    }
  }
  
  # Re-check after corrections
  if (is.data.frame(corrected_data)) {
    final_units <- unique(corrected_data[[mgmt_col]][!is.na(corrected_data[[mgmt_col]]) & corrected_data[[mgmt_col]] != ""])
  } else {
    final_units <- unique(corrected_data[!is.na(corrected_data) & corrected_data != ""])
  }
  
  final_missing <- setdiff(expected_units, final_units)
  final_unexpected <- setdiff(final_units, expected_units)
  
  # Create validation report
  validation_report <- list(
    analysis_name = analysis_name,
    expected_count = n_expected,
    found_count = length(final_units),
    missing_count = length(final_missing),
    unexpected_count = length(final_unexpected),
    missing_units = final_missing,
    unexpected_units = final_unexpected,
    corrections_made = corrections_made,
    all_units_present = length(final_missing) == 0,
    data_is_valid = length(final_missing) == 0 && length(final_unexpected) == 0
  )
  
  # Print validation results
  cat("\n=== MANAGEMENT UNIT VALIDATION:", analysis_name, "===\n")
  cat("Expected units:", n_expected, "| Found units:", length(final_units), "\n")
  
  if (length(corrections_made) > 0) {
    cat("Name corrections applied:\n")
    for (correction in corrections_made) {
      cat("  ", correction, "\n")
    }
  }
  
  if (length(final_missing) > 0) {
    cat("⚠️ MISSING UNITS (", length(final_missing), "):\n")
    for (unit in final_missing) {
      cat("  ❌", unit, "\n")
    }
  }
  
  if (length(final_unexpected) > 0) {
    cat("⚠️ UNEXPECTED UNITS (", length(final_unexpected), "):\n")
    for (unit in final_unexpected) {
      cat("  ❓", unit, "\n")
    }
  }
  
  if (validation_report$all_units_present) {
    cat("✅ All 21 management units are present!\n")
  } else {
    warning("Missing ", length(final_missing), " management units in ", analysis_name)
  }
  
  return(list(
    validation = validation_report,
    corrected_data = corrected_data,
    final_units = final_units
  ))
}

#' Apply watershed ordering to a data frame with validation
#' 
#' @param data Data frame containing management unit column
#' @param mgmt_col Character: name of management unit column
#' @param reverse_for_plots Logical: reverse order for coord_flip() plots
#' @param validate Logical: whether to run validation first
#' @param analysis_name Character: name of analysis for validation messages
#' @return Data frame with ordered factor
apply_watershed_order <- function(data, mgmt_col = "mgmt_river", reverse_for_plots = FALSE,
                                  validate = TRUE, analysis_name = "Analysis") {
  
  # Validate management units first if requested
  if (validate) {
    validation_result <- validate_management_units(data, mgmt_col, analysis_name)
    data <- validation_result$corrected_data
    
    # Stop if critical validation issues
    if (!validation_result$validation$all_units_present) {
      stop("Cannot proceed with ", analysis_name, " - missing management units. ",
           "Expected 21, found ", validation_result$validation$found_count)
    }
  }
  
  # Apply ordering
  standard_order <- get_watershed_order()
  if (reverse_for_plots) standard_order <- rev(standard_order)
  
  # Only include units that are actually in the data
  units_in_data <- unique(data[[mgmt_col]][!is.na(data[[mgmt_col]]) & data[[mgmt_col]] != ""])
  final_order <- standard_order[standard_order %in% units_in_data]
  
  # Apply factor ordering
  data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
  
  if (validate) {
    cat("\n=== FINAL WATERSHED ORDERING ===\n")
    cat("Direction:", ifelse(reverse_for_plots, "Downstream → Upstream (for plots)", "Upstream → Downstream"), "\n")
    for (i in 1:length(final_order)) {
      cat(sprintf("%2d. %s\n", i, final_order[i]))
    }
    cat("=================================\n")
  }
  
  return(data)
}

# =============================================================================
# QUARTILE DIVISION FUNCTIONS  
# =============================================================================

#' Divide natal data into DOY quartiles using fixed dates
#' 
#' @param natal_data Data frame with natal origins data
#' @return List with quartile subsets, breaks, and labels
divide_doy_quartiles <- function(natal_data) {
  
  # Fixed break dates (using 2024 for consistent DOY calculation)
  june_11 <- as.Date("2024-06-11")
  june_21 <- as.Date("2024-06-21")  
  july_01 <- as.Date("2024-07-01")
  
  june_11_doy <- as.numeric(format(june_11, "%j"))  # 163
  june_21_doy <- as.numeric(format(june_21, "%j"))  # 173
  july_01_doy <- as.numeric(format(july_01, "%j"))  # 183
  
  # Create quartile subsets using DOY
  quartile_subsets <- list(
    natal_data %>% filter(DOY <= june_11_doy),                    # Q1
    natal_data %>% filter(DOY > june_11_doy & DOY <= june_21_doy), # Q2
    natal_data %>% filter(DOY > june_21_doy & DOY <= july_01_doy), # Q3
    natal_data %>% filter(DOY > july_01_doy)                      # Q4
  )
  
  # Create labels
  subset_labels <- c(
    sprintf("Q1: Start-%d (to Jun 11)", june_11_doy),
    sprintf("Q2: %d-%d (Jun 12-21)", june_11_doy + 1, june_21_doy),
    sprintf("Q3: %d-%d (Jun 22-Jul 1)", june_21_doy + 1, july_01_doy),
    sprintf("Q4: %d-End (Jul 2+)", july_01_doy + 1)
  )
  
  return(list(
    subsets = quartile_subsets,
    breaks = c(min(natal_data$DOY, na.rm = TRUE), june_11_doy, june_21_doy, july_01_doy, max(natal_data$DOY, na.rm = TRUE)),
    labels = subset_labels
  ))
}

# =============================================================================
# BAYESIAN ASSIGNMENT FUNCTIONS
# =============================================================================

#' Calculate error values for Bayesian assignment
#' 
#' @param pid_isose Vector of isotope prediction standard errors
#' @param min_error Numeric: minimum error value to constrain predictions
#' @return Vector of calculated error values
calculate_error <- function(pid_isose, min_error) {
  
  # Constrain prediction errors above minimum
  pid_isose_mod <- ifelse(pid_isose < min_error, min_error, pid_isose)
  
  # Variance components
  within_site <- 0.0003133684 / 1.96
  analyt <- 0.00011 / 2
  
  # Combined error calculation
  error <- sqrt(pid_isose_mod^2 + within_site^2 + analyt^2)
  
  return(error)
}

#' Set up watershed-specific priors for assignment
#' 
#' @param edges SF object with stream edges
#' @param min_stream_order Numeric: minimum stream order
#' @param watershed Character: "Kusko" or "Yukon"
#' @param natal_data Data frame with natal data (unused for Kusko)
#' @return List of prior vectors
setup_watershed_priors <- function(edges, min_stream_order, watershed, natal_data = NULL) {
  
  StreamOrderPrior <- ifelse(edges$Str_Order >= min_stream_order, 1, 0)
  
  if (watershed == "Kusko") {
    pid_prior <- edges$UniPh2oNoE
    PresencePrior <- ifelse((edges$Str_Order %in% c(6, 7, 8)) & edges$SPAWNING_C == 0, 0, 1)
    
    return(list(
      pid_prior = pid_prior,
      StreamOrderPrior = StreamOrderPrior,
      PresencePrior = PresencePrior
    ))
    
  } else if (watershed == "Yukon") {
    # Simplified Yukon implementation
    pid_prior <- edges$PriorSl2
    PresencePrior <- ifelse((edges$Str_Order %in% c(7, 8, 9)) & edges$SPAWNING_C == 0, 0, 1)
    
    return(list(
      pid_prior = pid_prior,
      StreamOrderPrior = StreamOrderPrior,
      PresencePrior = PresencePrior
    ))
  }
  
  stop("Watershed must be 'Kusko' or 'Yukon'")
}

#' Perform Bayesian assignment for natal origins data
#' 
#' @param natal_data Data frame with natal origins data
#' @param edges SF object with stream edges  
#' @param watershed Character: "Kusko" or "Yukon"
#' @param priors List of prior values from setup_watershed_priors()
#' @param pid_iso Vector of isotope predictions for each edge
#' @param error Vector of error values for each edge
#' @param sensitivity_threshold Threshold for assignment filtering (default 0.7)
#' @return Matrix of assignments for each fish to each stream edge
perform_assignment <- function(natal_data, edges, watershed, priors, 
                               pid_iso, error, sensitivity_threshold = 0.7) {
  
  # Initialize assignment matrix (edges × fish)
  assignment_matrix <- matrix(NA, nrow = length(pid_iso), ncol = nrow(natal_data))
  
  # Process each fish individually
  for (i in 1:nrow(natal_data)) {
    iso_o <- as.numeric(natal_data$natal_iso[i])
    
    # Calculate assignment probabilities
    assign <- (1/sqrt(2*pi*error^2)) * exp(-1*(iso_o - pid_iso)^2/(2*error^2)) * 
      priors$pid_prior * priors$StreamOrderPrior * priors$PresencePrior
    
    # Normalize and threshold assignments
    assign_norm <- assign / sum(assign)
    assign_rescaled <- assign_norm / max(assign_norm)
    assign_rescaled[assign_rescaled < sensitivity_threshold] <- 0
    
    # Weight by CPUE
    assignment_matrix[,i] <- assign_rescaled * as.numeric(natal_data$COratio[i])
  }
  
  return(assignment_matrix)
}

# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================

#' Process management river data from assignment results with validation
#' 
#' @param edges SF object with stream edges
#' @param basin_assign_rescale Vector of assignment values
#' @param analysis_name Character: name of analysis for validation messages
#' @return Data frame with management unit summaries
process_mgmt_river_data <- function(edges, basin_assign_rescale, analysis_name = "Analysis") {
  
  if (!"mgmt_river" %in% colnames(edges)) {
    stop("mgmt_river column not found in edges data for ", analysis_name)
  }
  
  edges$basin_assign_rescale <- basin_assign_rescale
  
  # Filter to managed edges only
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  if (nrow(managed_edges) == 0) {
    stop("No edges with management unit assignments found in ", analysis_name)
  }
  
  # Summarize by management river
  summary_mgmt <- managed_edges %>%
    st_drop_geometry() %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(basin_assign_rescale, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    ) %>%
    mutate(production_proportion = total_production / sum(total_production, na.rm = TRUE))
  
  # Validate and ensure all 21 management units are present
  validation_result <- validate_management_units(summary_mgmt, "mgmt_river", analysis_name)
  summary_mgmt <- validation_result$corrected_data
  
  # Add missing units with zero values if needed
  if (!validation_result$validation$all_units_present) {
    missing_units <- validation_result$validation$missing_units
    
    missing_rows <- data.frame(
      mgmt_river = missing_units,
      total_production = 0,
      edge_count = 0,
      production_proportion = 0,
      stringsAsFactors = FALSE
    )
    
    summary_mgmt <- rbind(summary_mgmt, missing_rows)
    cat("✅ Added", length(missing_units), "missing units with zero production\n")
  }
  
  # Apply watershed ordering
  summary_mgmt <- apply_watershed_order(summary_mgmt, "mgmt_river", 
                                        reverse_for_plots = FALSE, 
                                        validate = FALSE,  # Already validated above
                                        analysis_name = analysis_name)
  
  return(summary_mgmt)
}

# =============================================================================
# PLOTTING UTILITY FUNCTIONS
# =============================================================================

#' Create output directory for analysis
#' 
#' @param analysis_name Character: name of analysis 
#' @param subdir Character: subdirectory ("figures" or "data")
#' @param base_dir Character: base output directory
#' @return Character path to directory
create_output_dir <- function(analysis_name, subdir = NULL, base_dir = "output") {
  if (is.null(subdir)) {
    dir_path <- file.path(base_dir, analysis_name)
  } else {
    dir_path <- file.path(base_dir, analysis_name, subdir)
  }
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  return(dir_path)
}

#' Save plot with consistent formatting
#' 
#' @param plot ggplot object or grob
#' @param filename Character: filename with extension
#' @param width Numeric: width in inches
#' @param height Numeric: height in inches
#' @param dpi Numeric: resolution
save_plot <- function(plot, filename, width = 12, height = 10, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white")
  cat("✅ Saved:", basename(filename), "\n")
}

#' Save data with consistent formatting
#' 
#' @param data Data frame to save
#' @param filename Character: filename with .csv extension
save_data <- function(data, filename) {
  write.csv(data, filename, row.names = FALSE)
  cat("✅ Saved:", basename(filename), "\n")
}