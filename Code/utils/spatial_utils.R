# spatial_utils.R
# Utility functions for spatial data processing in watershed analysis

library(sf)
library(dplyr)
library(here)
library(tidyr)

#' Process HUC data for watershed analysis
#' 
#' @param edges SF object containing stream edges
#' @param basin SF object containing basin boundary
#' @param Huc SF object containing HUC polygons
#' @param basin_assign_rescale Vector of assignment values for each edge
#' @param HUC Numeric HUC level (e.g., 8, 10)
#' @return SF object with HUC polygons and calculated metrics
process_huc_data <- function(edges, basin, Huc, basin_assign_rescale, HUC = 8) {
  huc_col <- paste0("HUC", HUC)
  name_col <- "Name"
  
  # remove the column HUC8 from edges 
  if (huc_col %in% colnames(edges)) {
    edges <- edges %>% select(-all_of(huc_col))
  }
  
  # Transform to consistent CRS
  edges <- st_transform(edges, st_crs(Huc))
  edges$basin_assign_rescale <- basin_assign_rescale
  basin <- st_transform(basin, st_crs(Huc))
  
  # Identify which HUCs intersect with the basin
  basin_buffer <- st_buffer(basin, dist = 0)
  hucs_in_basin <- Huc[st_intersects(Huc, basin_buffer, sparse = FALSE)[,1], ]
  
  # Calculate intersection areas
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
  
  # Calculate overlap percentage
  overlap_percentage <- intersection_areas %>%
    left_join(hucs_areas, by = huc_col) %>%
    mutate(pct_overlap = as.numeric(int_area / total_area))
  
  # Filter to include only HUCs with significant overlap
  significant_hucs <- overlap_percentage %>%
    filter(pct_overlap > 0.1) %>%
    pull(!!sym(huc_col))
  
  # Spatial join between edges and HUCs
  Combined_edges_HUC <- st_join(edges, Huc, join = st_intersects)
  edges$stream_length_m <- as.numeric(st_length(edges))
  
  # Summarize production by HUC polygon
  summary_huc <- Combined_edges_HUC %>%
    group_by(!!sym(huc_col)) %>%
    summarise(
      total_production = sum(basin_assign_rescale, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(production_proportion = total_production / sum(total_production, na.rm = TRUE))
  
  # Summarize stream length by HUC
  stream_length_by_huc <- edges %>%
    st_join(Huc, join = st_intersects) %>%
    st_drop_geometry() %>%
    group_by(!!sym(huc_col)) %>%
    summarise(total_stream_length = sum(stream_length_m, na.rm = TRUE))
  
  # Merge production and stream length data with HUC polygons
  final_result <- Huc %>%
    filter(!!sym(huc_col) %in% significant_hucs) %>%
    left_join(st_drop_geometry(summary_huc), by = huc_col) %>%
    left_join(stream_length_by_huc, by = huc_col)
  
  # Replace NA values with 0 using base R instead of replace_na
  final_result$total_production[is.na(final_result$total_production)] <- 0
  final_result$production_proportion[is.na(final_result$production_proportion)] <- 0
  final_result$total_stream_length[is.na(final_result$total_stream_length)] <- 0
  
  final_result <- final_result %>%
    mutate(
      production_per_meter = ifelse(total_stream_length > 0,
                                    total_production / total_stream_length,
                                    NA_real_),
      production_per_meter_norm = production_per_meter / max(production_per_meter, na.rm = TRUE)
    )
  
  return(final_result)
}

#' Process management unit data for watershed analysis
#' 
#' @param edges SF object containing stream edges with mgmt_river attribute
#' @param basin_assign_rescale Vector of assignment values for each edge
#' @param grand_total_production Total production across entire watershed for percentage calculations
#' @return Data frame with management unit summaries and calculated metrics
process_management_data <- function(edges, basin_assign_rescale, grand_total_production = NULL) {
  
  # Add the assignment values to edges
  edges$basin_assign_rescale <- basin_assign_rescale
  
  # Check if mgmt_river column exists
  if (!"mgmt_river" %in% colnames(edges)) {
    warning("mgmt_river column not found in edges data. Management analysis will be skipped.")
    return(NULL)
  }
  
  # Filter to only edges that have management unit assignments (not NA)
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  # Debug info
  total_edges <- nrow(edges)
  managed_count <- nrow(managed_edges)
  total_production_all <- sum(basin_assign_rescale, na.rm = TRUE)
  managed_production <- sum(managed_edges$basin_assign_rescale, na.rm = TRUE)
  
  message(paste("    Total edges:", total_edges, "| Managed edges:", managed_count, 
                "| Managed percentage:", round(100 * managed_count / total_edges, 1), "%"))
  message(paste("    Managed production as % of total:", 
                round(100 * managed_production / total_production_all, 1), "%"))
  
  if (managed_count == 0) {
    message("    No edges with management unit assignments found")
    return(NULL)
  }
  
  # Show unique management rivers
  unique_mgmt <- unique(managed_edges$mgmt_river)
  message(paste("    Management rivers found:", paste(unique_mgmt, collapse = ", ")))
  
  # Calculate stream length for each edge
  managed_edges$stream_length_m <- as.numeric(st_length(managed_edges))
  
  # Use provided grand total or calculate from current data
  if (is.null(grand_total_production)) {
    grand_total_production <- total_production_all
  }
  
  # Summarize production by management unit
  summary_mgmt <- managed_edges %>%
    st_drop_geometry() %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(basin_assign_rescale, na.rm = TRUE),
      total_stream_length = sum(stream_length_m, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    ) %>%
    # Calculate production proportion relative to ENTIRE RUN
    mutate(
      percent_of_total_run = (total_production / grand_total_production) * 100,
      production_per_meter = ifelse(total_stream_length > 0,
                                    total_production / total_stream_length,
                                    NA_real_)
    ) %>%
    # Calculate normalized production per meter
    mutate(
      production_per_meter_norm = ifelse(is.na(production_per_meter), 0,
                                         production_per_meter / max(production_per_meter, na.rm = TRUE))
    ) %>%
    # Add reference values
    mutate(
      grand_total_production = grand_total_production
    ) %>%
    # Handle any remaining NA values
    mutate(
      production_per_meter_norm = ifelse(is.na(production_per_meter_norm), 0, production_per_meter_norm)
    )
  
  return(summary_mgmt)
}

#' Load spatial data for a specific watershed
#'
#' @param watershed Character: "Kusko" or "Yukon"
#' @param HUC Numeric HUC level (e.g., 8, 10)
#' @param min_stream_order Minimum stream order to include
#' @return List containing edges, basin, and HUC data
load_spatial_data <- function(watershed, HUC = 8, min_stream_order = 3) {
  if (watershed == "Kusko") {
    # Use the ORIGINAL Kusko shapefile
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
  
  return(list(
    edges = edges,
    basin = basin,
    Huc = Huc
  ))
}

#' Load natal origins data
#'
#' @param year Character or numeric year
#' @param watershed Character: "Kusko" or "Yukon"
#' @return Data frame with natal origins data
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

#' Calculate error values for assignment
#'
#' @param pid_isose Vector of prediction error values
#' @param min_error Minimum error value to use
#' @return Vector of combined error values
calculate_error <- function(pid_isose, min_error) {
  # Constrain all pid_isose values above the minimum error value
  pid_isose_mod <- ifelse(pid_isose < min_error, min_error, pid_isose)
  
  # Calculate variance components
  within_site <- 0.0003133684 / 1.96
  analyt <- 0.00011 / 2
  
  # Calculate combined error
  error <- sqrt(pid_isose_mod^2 + within_site^2 + analyt^2)
  
  return(error)
}

#' Set up watershed-specific priors
#'
#' @param edges SF object with stream edges
#' @param min_stream_order Minimum stream order to include
#' @param watershed Character: "Kusko" or "Yukon"
#' @param natal_data Data frame with natal origins data (needed for Yukon)
#' @return List with priors and indices needed for assignments
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
    
    # Load Yukon-specific genetic groups
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

#' Divide data into quartiles based on fixed DOY dates (UPDATED)
#'
#' @param natal_data Data frame with natal origins data
#' @return List containing quartile subsets, breaks, and labels
divide_doy_quartiles <- function(natal_data) {
  # Fixed DOY breaks based on calendar dates
  # June 11 = DOY 162, then 10-day intervals
  # DOY 162 = June 11
  # DOY 172 = June 21 (162 + 10)
  # DOY 182 = July 1 (162 + 20) 
  # DOY 192 = July 11 (162 + 30)
  
  doy_breaks <- c(162, 172, 182, 192)
  
  # Create quartile subsets based on fixed dates
  quartile_subsets <- list()
  
  # Q1: Before June 11 (DOY < 162)
  quartile_subsets[[1]] <- natal_data %>% 
    filter(DOY < doy_breaks[1])
  
  # Q2: June 11-20 (DOY 162-171)  
  quartile_subsets[[2]] <- natal_data %>% 
    filter(DOY >= doy_breaks[1] & DOY < doy_breaks[2])
  
  # Q3: June 21-30 (DOY 172-181)
  quartile_subsets[[3]] <- natal_data %>% 
    filter(DOY >= doy_breaks[2] & DOY < doy_breaks[3])
  
  # Q4: July 1+ (DOY 182+)
  quartile_subsets[[4]] <- natal_data %>% 
    filter(DOY >= doy_breaks[3])
  
  # Create descriptive labels with actual date ranges found in data
  subset_labels <- character(4)
  
  for (i in 1:4) {
    if (nrow(quartile_subsets[[i]]) > 0) {
      actual_min <- min(quartile_subsets[[i]]$DOY, na.rm = TRUE)
      actual_max <- max(quartile_subsets[[i]]$DOY, na.rm = TRUE)
      
      # Convert DOY to readable dates for labels
      min_date <- format(as.Date(actual_min, origin = "2020-01-01"), "%b %d")
      max_date <- format(as.Date(actual_max, origin = "2020-01-01"), "%b %d")
      
      subset_labels[i] <- sprintf("DOY Q%d: %d-%d (%s-%s)", 
                                  i, actual_min, actual_max, min_date, max_date)
    } else {
      # Handle empty quartiles
      expected_ranges <- c("Before Jun 11", "Jun 11-20", "Jun 21-30", "Jul 1+")
      subset_labels[i] <- sprintf("DOY Q%d: %s (No data)", i, expected_ranges[i])
    }
  }
  
  return(list(
    subsets = quartile_subsets,
    breaks = doy_breaks,
    labels = subset_labels
  ))
}