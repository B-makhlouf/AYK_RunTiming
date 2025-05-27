################################################################################
# SPATIAL_UTILS.R - SPATIAL DATA PROCESSING UTILITIES
################################################################################
# PURPOSE: Core functions for loading and processing spatial data
# KEY CHANGE: HUC maps now show production PROPORTION, not scaled per stream length
################################################################################

library(sf)
library(dplyr)
library(here)
library(tidyr)

################################################################################
# HUC DATA PROCESSING - PRODUCTION PROPORTION (NOT PER STREAM LENGTH)
################################################################################

#' Process HUC data showing production proportion within each HUC
#' KEY: Shows proportion of production, NOT scaled to stream length
process_huc_data <- function(edges, basin, Huc, basin_assign_rescale, HUC = 8) {
  huc_col <- paste0("HUC", HUC)
  name_col <- "Name"
  
  # Remove existing HUC column from edges if present
  if (huc_col %in% colnames(edges)) {
    edges <- edges %>% select(-all_of(huc_col))
  }
  
  # Transform to consistent CRS and add assignment values
  edges <- st_transform(edges, st_crs(Huc))
  edges$basin_assign_rescale <- basin_assign_rescale
  basin <- st_transform(basin, st_crs(Huc))
  
  # Identify HUCs that intersect with basin
  basin_buffer <- st_buffer(basin, dist = 0)
  hucs_in_basin <- Huc[st_intersects(Huc, basin_buffer, sparse = FALSE)[,1], ]
  
  # Calculate intersection areas for filtering
  intersection_areas <- st_intersection(hucs_in_basin, basin_buffer) %>%
    mutate(area = st_area(.)) %>%
    st_drop_geometry() %>%
    group_by(!!sym(huc_col)) %>%
    summarize(int_area = sum(area))
  
  # Get original HUC areas
  hucs_areas <- hucs_in_basin %>%
    mutate(total_area = st_area(.)) %>%
    st_drop_geometry() %>%
    select(!!sym(huc_col), total_area)
  
  # Calculate overlap percentage and filter
  overlap_percentage <- intersection_areas %>%
    left_join(hucs_areas, by = huc_col) %>%
    mutate(pct_overlap = as.numeric(int_area / total_area))
  
  significant_hucs <- overlap_percentage %>%
    filter(pct_overlap > 0.1) %>%
    pull(!!sym(huc_col))
  
  # Spatial join between edges and HUCs
  Combined_edges_HUC <- st_join(edges, Huc, join = st_intersects)
  
  # CORE CALCULATION: Summarize production by HUC polygon (PROPORTION ONLY)
  summary_huc <- Combined_edges_HUC %>%
    group_by(!!sym(huc_col)) %>%
    summarise(
      total_production = sum(basin_assign_rescale, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(production_proportion = total_production / sum(total_production, na.rm = TRUE))
  
  # Merge with HUC polygons
  final_result <- Huc %>%
    filter(!!sym(huc_col) %in% significant_hucs) %>%
    left_join(st_drop_geometry(summary_huc), by = huc_col)
  
  # Replace NA values with 0
  final_result$total_production[is.na(final_result$total_production)] <- 0
  final_result$production_proportion[is.na(final_result$production_proportion)] <- 0
  
  return(final_result)
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
  
  return(summary_mgmt)
}

################################################################################
# DATA LOADING FUNCTIONS
################################################################################

#' Load spatial data for watershed analysis
load_spatial_data <- function(watershed, HUC = 8, min_stream_order = 3) {
  if (watershed == "Kusko") {
    edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC.shp", quiet = TRUE)
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
  
  Huc <- st_read(paste0("/Users/benjaminmakhlouf/Spatial Data/HUC", HUC, "_Trimmed.shp"), quiet = TRUE)
  Huc <- st_transform(Huc, st_crs(basin))
  
  return(list(edges = edges, basin = basin, Huc = Huc))
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

#' Divide natal data into DOY quartiles
divide_doy_quartiles <- function(natal_data) {
  # Get DOY range and create even breaks
  full_doy_range <- range(natal_data$DOY, na.rm = TRUE)
  doy_breaks <- seq(full_doy_range[1], full_doy_range[2], length.out = 5)
  
  # Create quartile subsets
  quartile_subsets <- list()
  for (i in 1:4) {
    quartile_subsets[[i]] <- natal_data %>% 
      filter(DOY >= doy_breaks[i] & DOY < doy_breaks[i+1])
  }
  
  # Create descriptive labels
  subset_labels <- sapply(1:4, function(i) {
    sprintf("DOY Q%d: %d-%d", i, ceiling(doy_breaks[i]), floor(doy_breaks[i+1]))
  })
  
  return(list(
    subsets = quartile_subsets,
    breaks = doy_breaks,
    labels = subset_labels
  ))
}