################################################################################
# SPATIAL_UTILS.R - SPATIAL DATA PROCESSING UTILITIES
################################################################################
# PURPOSE: Core functions for loading and processing spatial data
# NOTE: HUC processing removed - focus on management river analysis only
################################################################################

library(sf)
library(dplyr)
library(here)
library(tidyr)

################################################################################
# WATERSHED ORDERING FUNCTIONS
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
  
  # Add any units from data that aren't in our standard list (shouldn't happen but safety check)
  missing_units <- setdiff(units_in_data, standard_order)
  if (length(missing_units) > 0) {
    warning("Found management units not in standard order: ", paste(missing_units, collapse = ", "))
    final_order <- c(final_order, missing_units)
  }
  
  # Apply factor ordering
  data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
  
  return(data)
}

################################################################################
# MANAGEMENT RIVER DATA PROCESSING - PRODUCTION PROPORTION
################################################################################

#' Process management river data showing production proportion within each management unit
#' Similar to HUC processing but groups by mgmt_river attribute
process_mgmt_river_data <- function(edges, basin_assign_rescale) {
  
  # Check if management data exists
  if (!"mgmt_river" %in% colnames(edges)) {
    warning("mgmt_river column not found in edges data. Management analysis will be skipped.")
    return(NULL)
  }
  
  # Add assignment values to edges
  edges$basin_assign_rescale <- basin_assign_rescale
  
  # Filter to only edges that have management unit assignments (not NA or empty)
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  if (nrow(managed_edges) == 0) {
    message("No edges with management unit assignments found")
    return(NULL)
  }
  
  message(paste("Found", nrow(managed_edges), "managed edges out of", nrow(edges), "total edges"))
  
  # CORE CALCULATION: Summarize production by management river (PROPORTION ONLY)
  summary_mgmt <- managed_edges %>%
    st_drop_geometry() %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(basin_assign_rescale, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    ) %>%
    mutate(production_proportion = total_production / sum(total_production, na.rm = TRUE))
  
  # Apply watershed ordering for consistency
  summary_mgmt <- apply_watershed_order(summary_mgmt, "mgmt_river", reverse_for_plots = FALSE)
  
  return(summary_mgmt)
}

################################################################################
# MANAGEMENT UNIT DATA PROCESSING
################################################################################

#' Process management unit data for watershed analysis
process_management_data <- function(edges, basin_assign_rescale, grand_total_production = NULL) {
  
  edges$basin_assign_rescale <- basin_assign_rescale
  
  # Check if management data exists
  if (!"mgmt_river" %in% colnames(edges)) {
    warning("mgmt_river column not found in edges data. Management analysis will be skipped.")
    return(NULL)
  }
  
  # Filter to managed edges only
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  if (nrow(managed_edges) == 0) {
    message("No edges with management unit assignments found")
    return(NULL)
  }
  
  # Calculate stream length and use provided total or calculate
  managed_edges$stream_length_m <- as.numeric(st_length(managed_edges))
  
  if (is.null(grand_total_production)) {
    grand_total_production <- sum(basin_assign_rescale, na.rm = TRUE)
  }
  
  # Summarize by management unit
  summary_mgmt <- managed_edges %>%
    st_drop_geometry() %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(basin_assign_rescale, na.rm = TRUE),
      total_stream_length = sum(stream_length_m, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      percent_of_total_run = (total_production / grand_total_production) * 100,
      grand_total_production = grand_total_production
    )
  
  # Apply watershed ordering for consistency
  summary_mgmt <- apply_watershed_order(summary_mgmt, "mgmt_river", reverse_for_plots = FALSE)
  
  return(summary_mgmt)
}

################################################################################
# DATA LOADING FUNCTIONS
################################################################################

#' Load spatial data for watershed analysis
load_spatial_data <- function(watershed, HUC = 8, min_stream_order = 3) {
  if (watershed == "Kusko") {
    edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp", quiet = TRUE)
    basin <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp", quiet = TRUE)
  } else if (watershed == "Yukon") {
    edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/USGS Added/YukonUSGS.shp", quiet = TRUE)
    basin <- st_read("/Users/benjaminmakhlouf/Spatial Data/Basin Map Necessary Shapefiles/Yuk_Mrg_final_alb.shp", quiet = TRUE)
  } else {
    stop("Watershed must be either 'Kusko' or 'Yukon'")
  }
  
  # Transform and filter data
  edges <- st_transform(edges, st_crs(basin))
  edges <- edges[edges$Str_Order >= min_stream_order,]
  
  # Note: HUC loading removed since we're not using HUC maps anymore
  # Return without HUC data
  return(list(edges = edges, basin = basin))
}

#' Load natal origins data for specific year and watershed
load_natal_data <- function(year, watershed) {
  file_path <- paste0(
    "/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/Data/",
    "Natal Origin Analysis Data/03_Natal Origins Genetics CPUE/",
    year, "_", watershed, "_Natal_Origins_Genetics_CPUE.csv"
  )
  
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
  
  return(clean_data)
}

################################################################################
# ASSIGNMENT CALCULATION UTILITIES
################################################################################

#' Calculate error values for Bayesian assignment
calculate_error <- function(pid_isose, min_error) {
  # Constrain prediction errors above minimum
  pid_isose_mod <- ifelse(pid_isose < min_error, min_error, pid_isose)
  
  # Calculate variance components
  within_site <- 0.0003133684 / 1.96
  analyt <- 0.00011 / 2
  
  # Combined error calculation
  error <- sqrt(pid_isose_mod^2 + within_site^2 + analyt^2)
  
  return(error)
}

#' Set up watershed-specific priors for assignment
setup_watershed_priors <- function(edges, min_stream_order, watershed, natal_data = NULL) {
  StreamOrderPrior <- ifelse(edges$Str_Order >= min_stream_order, 1, 0)
  
  if (watershed == "Kusko") {
    pid_prior <- edges$UniPh2oNoE
    PresencePrior <- ifelse((edges$Str_Order %in% c(6, 7, 8)) & edges$SPAWNING_C == 0, 0, 1)
    NewHabitatPrior <- ifelse(edges$Spawner_IP == 0, 0, 1)
    
    return(list(
      pid_prior = pid_prior,
      StreamOrderPrior = StreamOrderPrior,
      PresencePrior = PresencePrior,
      NewHabitatPrior = NewHabitatPrior
    ))
    
  } else if (watershed == "Yukon") {
    pid_prior <- edges$PriorSl2
    
    # Load Yukon genetic groups
    ly.gen <- st_read(here("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Yukon/For_Sean/edges_LYGen.shp"), quiet = TRUE)
    ly.gen_reachid <- ly.gen$reachid
    my.gen <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Yukon/For_Sean/edges_MYGen.shp", quiet = TRUE)
    my.gen_reachid <- my.gen$reachid
    uy.gen <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Yukon/For_Sean/edges_UYGen.shp", quiet = TRUE)
    uy.gen_reachid <- uy.gen$reachid
    
    edges$GenLMU <- 0
    edges$GenLMU[edges$reachid %in% ly.gen_reachid] <- "lower"
    edges$GenLMU[edges$reachid %in% my.gen_reachid] <- "middle"
    edges$GenLMU[edges$reachid %in% uy.gen_reachid] <- "upper"
    
    LYsites <- which(edges$GenLMU == "lower")
    MYsites <- which(edges$GenLMU == "middle")
    UYsites <- which(edges$GenLMU == "upper")
    
    PresencePrior <- ifelse((edges$Str_Order %in% c(7, 8, 9)) & edges$SPAWNING_C == 0, 0, 1)
    NewHabitatPrior <- ifelse(edges$Spawner_IP == 0, 0, 1)
    
    return(list(
      pid_prior = pid_prior,
      StreamOrderPrior = StreamOrderPrior,
      PresencePrior = PresencePrior,
      NewHabitatPrior = NewHabitatPrior,
      LYsites = LYsites,
      MYsites = MYsites,
      UYsites = UYsites
    ))
  }
}

################################################################################
# DOY QUARTILE DIVISION
################################################################################

#' Divide natal data into custom DOY quartiles
divide_doy_quartiles <- function(natal_data) {
  # Handle different date formats
  if ("Date" %in% colnames(natal_data)) {
    # Try to convert Date column to proper Date format
    if (!inherits(natal_data$Date, "Date")) {
      natal_data$Date <- as.Date(natal_data$Date)
    }
    year <- format(natal_data$Date[1], "%Y")
  } else {
    # Fall back to extracting year from DOY (assuming current approach)
    # You'll need to specify which year this data represents
    warning("No Date column found, using DOY. Please specify year manually.")
    year <- "2020"  # Default - change this to match your data year
  }
  
  # Define fixed break dates for this year
  june_11 <- as.Date(paste0(year, "-06-11"))
  june_21 <- as.Date(paste0(year, "-06-21"))  
  july_01 <- as.Date(paste0(year, "-07-01"))
  
  # Create quartile subsets using actual dates
  if ("Date" %in% colnames(natal_data)) {
    quartile_subsets <- list(
      natal_data %>% filter(Date <= june_11),           # Q1: Start to Jun 11
      natal_data %>% filter(Date > june_11 & Date <= june_21), # Q2: Jun 12-21
      natal_data %>% filter(Date > june_21 & Date <= july_01), # Q3: Jun 22-Jul 1  
      natal_data %>% filter(Date > july_01)             # Q4: Jul 2-End
    )
  } else {
    # Fall back to DOY-based filtering
    june_11_doy <- as.numeric(format(june_11, "%j"))
    june_21_doy <- as.numeric(format(june_21, "%j"))
    july_01_doy <- as.numeric(format(july_01, "%j"))
    
    quartile_subsets <- list(
      natal_data %>% filter(DOY <= june_11_doy),           
      natal_data %>% filter(DOY > june_11_doy & DOY <= june_21_doy), 
      natal_data %>% filter(DOY > june_21_doy & DOY <= july_01_doy),  
      natal_data %>% filter(DOY > july_01_doy)             
    )
  }
  
  # Convert to DOY for labels
  june_11_doy <- as.numeric(format(june_11, "%j"))
  june_21_doy <- as.numeric(format(june_21, "%j"))
  july_01_doy <- as.numeric(format(july_01, "%j"))
  
  # Create labels
  subset_labels <- c(
    sprintf("DOY Q1: Start-%d (to Jun 11)", june_11_doy),
    sprintf("DOY Q2: %d-%d (Jun 12-21)", june_11_doy + 1, june_21_doy),
    sprintf("DOY Q3: %d-%d (Jun 22-Jul 1)", june_21_doy + 1, july_01_doy),
    sprintf("DOY Q4: %d-End (Jul 2+)", july_01_doy + 1)
  )
  
  # DOY breaks for compatibility
  doy_breaks <- c(
    min(natal_data$DOY, na.rm = TRUE),
    june_11_doy,
    june_21_doy, 
    july_01_doy,
    max(natal_data$DOY, na.rm = TRUE)
  )
  
  return(list(
    subsets = quartile_subsets,
    breaks = doy_breaks,
    labels = subset_labels
  ))
}