################################################################################
# SPATIAL_UTILS.R - COMPLETE SPATIAL UTILITIES WITH ALL MISSING FUNCTIONS
################################################################################
# PURPOSE: Contains all spatial data loading and utility functions needed for
#          the salmon run timing analysis
################################################################################

library(sf)
library(dplyr)

################################################################################
# CORRECTED WATERSHED ORDERING FUNCTIONS
################################################################################

#' Get standardized watershed ordering for all plots and analyses
#' 
#' @return Character vector of management unit names ordered from upstream to downstream
get_watershed_order <- function() {
  # Exact order from your document - upstream to downstream
  watershed_order <- c(
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
  
  return(watershed_order)
}

#' Apply watershed ordering to a data frame
#' 
#' @param data Data frame containing management unit column
#' @param mgmt_col Name of the management unit column (default: "mgmt_river")
#' @param reverse_for_plots Logical, whether to reverse order for coord_flip() plots (default: FALSE)
#' @return Data frame with management unit column converted to ordered factor
apply_watershed_order <- function(data, mgmt_col = "mgmt_river", reverse_for_plots = FALSE) {
  
  # Get the standard ordering
  standard_order <- get_watershed_order()
  
  # Reverse if needed for coord_flip plots (so upstream appears at top)
  if (reverse_for_plots) {
    standard_order <- rev(standard_order)
  }
  
  # Filter to only units present in the data
  units_in_data <- unique(data[[mgmt_col]])
  final_order <- standard_order[standard_order %in% units_in_data]
  
  # Handle name variations/mismatches
  name_mapping <- c(
    "E. Fork Kuskokwim" = "E. Fork Kuskokwim River"
  )
  
  # Apply any name mappings if needed
  for (old_name in names(name_mapping)) {
    if (old_name %in% units_in_data && !(name_mapping[old_name] %in% units_in_data)) {
      data[[mgmt_col]][data[[mgmt_col]] == old_name] <- name_mapping[old_name]
      units_in_data <- unique(data[[mgmt_col]])
      final_order <- standard_order[standard_order %in% units_in_data]
    }
  }
  
  # Add any units from data that aren't in our standard list
  missing_units <- setdiff(units_in_data, standard_order)
  if (length(missing_units) > 0) {
    warning("Found management units not in standard order: ", paste(missing_units, collapse = ", "))
    message("Missing units (will be added at end):")
    for(unit in missing_units) message("  - ", unit)
    final_order <- c(final_order, missing_units)
  }
  
  # Apply factor ordering
  data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
  
  return(data)
}

################################################################################
# SPATIAL DATA LOADING FUNCTIONS
################################################################################

#' Load spatial data for a given watershed
#' @param watershed Character: "Kusko" or "Yukon" 
#' @param huc_level Numeric: HUC level (typically 8)
#' @param min_stream_order Numeric: minimum stream order to include
#' @return List containing edges, basin, and Huc spatial data
load_spatial_data <- function(watershed, huc_level = 8, min_stream_order = 3) {
  
  cat("Loading spatial data for", watershed, "watershed...\n")
  
  if (watershed == "Kusko") {
    # Kuskokwim watershed file paths
    edges_path <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
    basin_path <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
    huc_path <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC.shp"
    
    if (!file.exists(edges_path)) {
      stop("Kusko edges file not found at: ", edges_path)
    }
    if (!file.exists(basin_path)) {
      stop("Kusko basin file not found at: ", basin_path)
    }
    
  } else if (watershed == "Yukon") {
    # Yukon watershed file paths - adjust these to your actual paths
    edges_path <- "/Users/benjaminmakhlouf/Spatial Data/Yukon_edges.shp"
    basin_path <- "/Users/benjaminmakhlouf/Spatial Data/Yukon_basin.shp"
    huc_path <- "/Users/benjaminmakhlouf/Spatial Data/Yukon_HUC.shp"
    
    if (!file.exists(edges_path)) {
      stop("Yukon edges file not found at: ", edges_path, "\nPlease update the path in load_spatial_data function")
    }
    if (!file.exists(basin_path)) {
      stop("Yukon basin file not found at: ", basin_path, "\nPlease update the path in load_spatial_data function")
    }
    
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon', got: ", watershed)
  }
  
  # Load the spatial data
  tryCatch({
    edges <- st_read(edges_path, quiet = TRUE)
    basin <- st_read(basin_path, quiet = TRUE)
    
    # Load HUC data if it exists
    if (file.exists(huc_path)) {
      huc_data <- st_read(huc_path, quiet = TRUE)
    } else {
      huc_data <- NULL
      warning("HUC file not found at: ", huc_path)
    }
    
    # Filter by stream order if specified
    if (!is.null(min_stream_order) && "Str_Order" %in% colnames(edges)) {
      original_count <- nrow(edges)
      edges <- edges %>% filter(Str_Order >= min_stream_order)
      cat("Filtered edges from", original_count, "to", nrow(edges), "features (stream order >=", min_stream_order, ")\n")
    }
    
    cat("✓ Successfully loaded spatial data for", watershed, "\n")
    cat("  - Edges:", nrow(edges), "features\n")
    cat("  - Basin: 1 polygon\n")
    if (!is.null(huc_data)) {
      cat("  - HUC data:", nrow(huc_data), "features\n")
    }
    
    return(list(
      edges = edges,
      basin = basin,
      Huc = huc_data
    ))
    
  }, error = function(e) {
    stop("Error loading spatial data for ", watershed, ": ", e$message)
  })
}

#' Load natal data for a specific year and watershed
#' @param year Numeric: year to load data for
#' @param watershed Character: "Kusko" or "Yukon"
#' @return Data frame with natal data
load_natal_data <- function(year, watershed) {
  
  cat("Loading natal data for", watershed, year, "...\n")
  
  if (watershed == "Kusko") {
    # Adjust this path to your actual natal data location
    data_path <- paste0("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/", watershed, "_", year, "_natal_data.csv")
    
    # Alternative path pattern if the above doesn't work
    if (!file.exists(data_path)) {
      data_path <- paste0("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Natal_Data/", watershed, "_", year, ".csv")
    }
    
    # Another alternative
    if (!file.exists(data_path)) {
      data_path <- paste0("/Users/benjaminmakhlouf/Data/", watershed, "/", year, "_natal_origins.csv")
    }
    
  } else if (watershed == "Yukon") {
    # Yukon data path - adjust as needed
    data_path <- paste0("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/", watershed, "_", year, "_natal_data.csv")
    
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon', got: ", watershed)
  }
  
  if (!file.exists(data_path)) {
    stop("Natal data file not found at: ", data_path, "\nPlease update the path in load_natal_data function")
  }
  
  tryCatch({
    natal_data <- read.csv(data_path)
    
    # Validate required columns
    required_cols <- c("DOY", "natal_iso", "COratio", "dailyCPUEprop")
    missing_cols <- setdiff(required_cols, colnames(natal_data))
    
    if (length(missing_cols) > 0) {
      stop("Missing required columns in natal data: ", paste(missing_cols, collapse = ", "))
    }
    
    cat("✓ Successfully loaded natal data:", nrow(natal_data), "observations\n")
    return(natal_data)
    
  }, error = function(e) {
    stop("Error loading natal data for ", watershed, " ", year, ": ", e$message)
  })
}

#' Divide natal data into DOY quartiles
#' @param natal_data Data frame with natal data containing DOY column
#' @return List with quartile subsets and labels
divide_doy_quartiles <- function(natal_data) {
  
  if (!"DOY" %in% colnames(natal_data)) {
    stop("DOY column not found in natal data")
  }
  
  # Calculate quartile breakpoints
  doy_quartiles <- quantile(natal_data$DOY, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  
  # Create quartile subsets
  q1_data <- natal_data %>% filter(DOY >= doy_quartiles[1] & DOY <= doy_quartiles[2])
  q2_data <- natal_data %>% filter(DOY > doy_quartiles[2] & DOY <= doy_quartiles[3])
  q3_data <- natal_data %>% filter(DOY > doy_quartiles[3] & DOY <= doy_quartiles[4])
  q4_data <- natal_data %>% filter(DOY > doy_quartiles[4] & DOY <= doy_quartiles[5])
  
  # Create labels
  labels <- c(
    paste0("Q1: DOY ", round(doy_quartiles[1]), "-", round(doy_quartiles[2])),
    paste0("Q2: DOY ", round(doy_quartiles[2]), "-", round(doy_quartiles[3])),
    paste0("Q3: DOY ", round(doy_quartiles[3]), "-", round(doy_quartiles[4])),
    paste0("Q4: DOY ", round(doy_quartiles[4]), "-", round(doy_quartiles[5]))
  )
  
  cat("Divided natal data into quartiles:\n")
  cat("  Q1:", nrow(q1_data), "observations\n")
  cat("  Q2:", nrow(q2_data), "observations\n")
  cat("  Q3:", nrow(q3_data), "observations\n")
  cat("  Q4:", nrow(q4_data), "observations\n")
  
  return(list(
    subsets = list(q1_data, q2_data, q3_data, q4_data),
    labels = labels,
    quartiles = doy_quartiles
  ))
}

#' Calculate error values for assignment
#' @param pid_isose Vector of isotope standard errors
#' @param min_error Minimum error value to use
#' @return Vector of error values
calculate_error <- function(pid_isose, min_error = 0.0006) {
  
  if (is.null(pid_isose)) {
    stop("pid_isose cannot be null")
  }
  
  # Use the larger of pid_isose or min_error
  error <- pmax(pid_isose, min_error, na.rm = TRUE)
  
  # Handle NA values
  error[is.na(error)] <- min_error
  
  return(error)
}

#' Setup watershed priors for assignment
#' @param edges Spatial edges data
#' @param min_stream_order Minimum stream order
#' @param watershed Watershed name
#' @param natal_data Natal data (for additional prior calculations if needed)
#' @return List of prior vectors
setup_watershed_priors <- function(edges, min_stream_order, watershed, natal_data = NULL) {
  
  cat("Setting up watershed priors for", watershed, "...\n")
  
  # Stream order prior
  if ("Str_Order" %in% colnames(edges)) {
    StreamOrderPrior <- ifelse(edges$Str_Order >= min_stream_order, 1, 0)
  } else {
    warning("Str_Order column not found in edges, using all streams")
    StreamOrderPrior <- rep(1, nrow(edges))
  }
  
  # Basic watershed prior (can be refined based on specific watershed knowledge)
  if ("iso_pred" %in% colnames(edges)) {
    pid_prior <- ifelse(!is.na(edges$iso_pred), 1, 0)
  } else {
    warning("iso_pred column not found in edges, using uniform prior")
    pid_prior <- rep(1, nrow(edges))
  }
  
  # Presence prior (streams that are present/accessible)
  PresencePrior <- rep(1, nrow(edges))  # Assume all streams are accessible unless specified otherwise
  
  # Watershed-specific priors
  if (watershed == "Yukon") {
    # For Yukon, we might have genetic regions
    # These would need to be defined based on your specific data
    # For now, using placeholder logic
    
    total_streams <- nrow(edges)
    LYsites <- 1:round(total_streams * 0.3)  # Lower Yukon sites (adjust as needed)
    MYsites <- round(total_streams * 0.3 + 1):round(total_streams * 0.7)  # Middle Yukon sites
    UYsites <- round(total_streams * 0.7 + 1):total_streams  # Upper Yukon sites
    
    cat("  - Lower Yukon sites:", length(LYsites), "\n")
    cat("  - Middle Yukon sites:", length(MYsites), "\n") 
    cat("  - Upper Yukon sites:", length(UYsites), "\n")
    
    return(list(
      StreamOrderPrior = StreamOrderPrior,
      pid_prior = pid_prior,
      PresencePrior = PresencePrior,
      LYsites = LYsites,
      MYsites = MYsites,
      UYsites = UYsites
    ))
    
  } else {
    # For Kuskokwim and other watersheds
    cat("  - Stream order prior: ", sum(StreamOrderPrior), "/", length(StreamOrderPrior), "streams included\n")
    cat("  - Watershed prior: ", sum(pid_prior), "/", length(pid_prior), "streams with data\n")
    
    return(list(
      StreamOrderPrior = StreamOrderPrior,
      pid_prior = pid_prior,
      PresencePrior = PresencePrior
    ))
  }
}

#' Process management river data from stream assignments
#' @param edges Spatial edges data with mgmt_river column
#' @param stream_assignments Vector of assignments for each stream
#' @return Data frame with management river summaries
process_mgmt_river_data <- function(edges, stream_assignments) {
  
  if (!"mgmt_river" %in% colnames(edges)) {
    warning("mgmt_river column not found in edges data")
    return(NULL)
  }
  
  if (length(stream_assignments) != nrow(edges)) {
    stop("Length of stream_assignments (", length(stream_assignments), 
         ") does not match number of edges (", nrow(edges), ")")
  }
  
  # Create data frame combining edges info with assignments
  mgmt_data <- data.frame(
    mgmt_river = edges$mgmt_river,
    assignment_value = stream_assignments,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(mgmt_river) & mgmt_river != "") %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(assignment_value, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    )
  
  # Calculate production proportion (within this dataset)
  total_all_production <- sum(mgmt_data$total_production, na.rm = TRUE)
  
  if (total_all_production > 0) {
    mgmt_data$production_proportion <- mgmt_data$total_production / total_all_production
  } else {
    mgmt_data$production_proportion <- 0
  }
  
  # Sort by production (highest first)
  mgmt_data <- mgmt_data %>%
    arrange(desc(total_production))
  
  cat("Processed management river data for", nrow(mgmt_data), "management units\n")
  
  return(mgmt_data)
}

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Get watershed-specific parameters
#' @param watershed Character: "Kusko" or "Yukon"
#' @return List of watershed-specific parameters
get_watershed_params <- function(watershed) {
  if (watershed == "Kusko") {
    list(
      sensitivity_threshold = 0.7, 
      min_error = 0.0006, 
      min_stream_order = 3
    )
  } else if (watershed == "Yukon") {
    list(
      sensitivity_threshold = 0.7, 
      min_error = 0.003, 
      min_stream_order = 5
    )
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon'")
  }
}

#' Check if required data files exist for a watershed and year
#' @param watershed Character: watershed name
#' @param year Numeric: year
#' @return Logical: TRUE if all files exist
check_data_availability <- function(watershed, year) {
  
  # Check spatial data
  if (watershed == "Kusko") {
    spatial_files <- c(
      "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp",
      "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
    )
  } else if (watershed == "Yukon") {
    spatial_files <- c(
      "/Users/benjaminmakhlouf/Spatial Data/Yukon_edges.shp",
      "/Users/benjaminmakhlouf/Spatial Data/Yukon_basin.shp"
    )
  } else {
    return(FALSE)
  }
  
  spatial_exist <- all(sapply(spatial_files, file.exists))
  
  # Check natal data (simplified - you may need to adjust paths)
  natal_path <- paste0("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/", watershed, "_", year, "_natal_data.csv")
  natal_exist <- file.exists(natal_path)
  
  if (!spatial_exist) {
    cat("Missing spatial files for", watershed, "\n")
  }
  if (!natal_exist) {
    cat("Missing natal data for", watershed, year, "\n")
  }
  
  return(spatial_exist && natal_exist)
}