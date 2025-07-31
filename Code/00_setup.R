################################################################################
# 00_SETUP.R - SALMON RUN TIMING ANALYSIS SETUP
################################################################################
# All parameters, paths, and essential functions for the analysis suite
# Load this first before running any other scripts
################################################################################

# Load all required libraries
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(tidyr)
  library(scales)
  library(readr)
  library(here)
  library(glue)
  library(grid)
  library(gridExtra)
  library(viridis)
})

cat("=== SALMON RUN TIMING ANALYSIS SETUP ===\n")

################################################################################
# GLOBAL PARAMETERS
################################################################################

# Watershed parameters
WATERSHED_PARAMS <- list(
  "Kusko" = list(
    sensitivity_threshold = 0.7,
    min_error = 0.0006,
    min_stream_order = 3
  ),
  "Yukon" = list(
    sensitivity_threshold = 0.7,
    min_error = 0.003,
    min_stream_order = 5
  )
)

# Analysis configuration
CONFIG <- list(
  years = c(2017, 2018, 2019, 2020, 2021, 2022),
  watersheds = c("Kusko"),
  quartiles = c("Q1", "Q2", "Q3", "Q4"),
  quartile_colors = c("Q1" = "#E31A1C", "Q2" = "#1F78B4", 
                      "Q3" = "#33A02C", "Q4" = "#FF7F00"),
  closure_dates = list(start = "06-01", end = "06-11"),
  cumulative_interval_days = 3
)

# File paths - UPDATE THE BASE_DIR FOR YOUR SYSTEM
BASE_DIR <- "/Users/benjaminmakhlouf"

PATHS <- list(
  # Spatial data
  kusko_edges = file.path(BASE_DIR, "Spatial Data/KuskoUSGS_HUC_joined.shp"),
  kusko_basin = file.path(BASE_DIR, "Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"),
  yukon_edges = file.path(BASE_DIR, "Spatial Data/USGS Added/YukonUSGS.shp"),
  yukon_basin = file.path(BASE_DIR, "Spatial Data/Basin Map Necessary Shapefiles/Yuk_Mrg_final_alb.shp"),
  
  # Natal data directory
  natal_data_dir = file.path(BASE_DIR, 
                             "Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/Data",
                             "Natal Origin Analysis Data/03_Natal Origins Genetics CPUE"),
  
  # Output directories
  output_dir = "Analysis_Results",
  maps_dir = "Basin_Maps", 
  figures_dir = "Figures"
)

# Management unit ordering (upstream to downstream)
# Management unit ordering (upstream to downstream) - CORRECTED
WATERSHED_ORDER <- c(
  "N. Fork Kusko", 
  "E. Fork Kuskokwim",        # East Fork in proper position
  "S. Fork Kusko",
  "Upper Kusko Main", 
  "Big River", 
  "Takotna and Nixon Fork",
  "Tatlawiksuk", 
  "Swift", 
  "Stony", 
  "Holitna",                  # SEPARATE unit
  "Hoholitna",                # SEPARATE unit (this was missing!)
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

################################################################################
# ESSENTIAL DATA LOADING FUNCTIONS
################################################################################

#' Load spatial data for a watershed
load_spatial_data <- function(watershed) {
  params <- WATERSHED_PARAMS[[watershed]]
  
  if (watershed == "Kusko") {
    edges <- st_read(PATHS$kusko_edges, quiet = TRUE)
    basin <- st_read(PATHS$kusko_basin, quiet = TRUE)
  } else if (watershed == "Yukon") {
    edges <- st_read(PATHS$yukon_edges, quiet = TRUE)
    basin <- st_read(PATHS$yukon_basin, quiet = TRUE)
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon'")
  }
  
  # Transform CRS and filter by stream order
  edges <- st_transform(edges, st_crs(basin))
  edges <- edges[edges$Str_Order >= params$min_stream_order, ]
  
  return(list(edges = edges, basin = basin))
}

#' Load natal origins data for a specific year and watershed
load_natal_data <- function(year, watershed) {
  file_path <- file.path(PATHS$natal_data_dir, 
                         paste0(year, "_", watershed, "_Natal_Origins_Genetics_CPUE.csv"))
  
  if (!file.exists(file_path)) {
    stop("Natal data file not found: ", file_path)
  }
  
  natal_data <- read_csv(file_path, show_col_types = FALSE)
  
  # Clean data based on watershed requirements
  if (watershed == "Yukon") {
    clean_data <- natal_data %>%
      filter(!is.na(Lower), !is.na(natal_iso), !is.na(dailyCPUEprop))
  } else {
    clean_data <- natal_data %>%
      filter(!is.na(natal_iso), !is.na(dailyCPUEprop))
  }
  
  return(clean_data)
}

#' Apply consistent watershed ordering to data
apply_watershed_order <- function(data, mgmt_col = "mgmt_river", reverse_for_plots = FALSE) {
  order_to_use <- if (reverse_for_plots) rev(WATERSHED_ORDER) else WATERSHED_ORDER
  
  units_in_data <- unique(data[[mgmt_col]])
  final_order <- order_to_use[order_to_use %in% units_in_data]
  
  # Add any units not in standard order
  missing_units <- setdiff(units_in_data, WATERSHED_ORDER)
  if (length(missing_units) > 0) {
    warning("Management units not in standard order: ", paste(missing_units, collapse = ", "))
    final_order <- c(final_order, missing_units)
  }
  
  data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
  return(data)
}

################################################################################
# CORE CALCULATION FUNCTIONS
################################################################################

#' Calculate error values for Bayesian assignment
calculate_error <- function(pid_isose, min_error) {
  pid_isose_mod <- ifelse(pid_isose < min_error, min_error, pid_isose)
  within_site <- 0.0003133684 / 1.96
  analyt <- 0.00011 / 2
  error <- sqrt(pid_isose_mod^2 + within_site^2 + analyt^2)
  return(error)
}

#' Set up watershed-specific priors
setup_priors <- function(edges, watershed, natal_data = NULL) {
  params <- WATERSHED_PARAMS[[watershed]]
  StreamOrderPrior <- ifelse(edges$Str_Order >= params$min_stream_order, 1, 0)
  
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
    PresencePrior <- ifelse((edges$Str_Order %in% c(7, 8, 9)) & edges$SPAWNING_C == 0, 0, 1)
    NewHabitatPrior <- ifelse(edges$Spawner_IP == 0, 0, 1)
    
    # Load Yukon genetic groups (simplified for this example)
    LYsites <- which(edges$GenLMU == "lower")
    MYsites <- which(edges$GenLMU == "middle") 
    UYsites <- which(edges$GenLMU == "upper")
    
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

#' Perform Bayesian assignment
perform_assignment <- function(natal_data, edges, watershed, priors, pid_iso, error, sensitivity_threshold) {
  assignment_matrix <- matrix(NA, nrow = length(pid_iso), ncol = nrow(natal_data))
  
  for (i in 1:nrow(natal_data)) {
    iso_o <- as.numeric(natal_data$natal_iso[i])
    
    if (watershed == "Kusko") {
      assign <- (1/sqrt(2*pi*error^2)) * exp(-1*(iso_o - pid_iso)^2/(2*error^2)) * 
        priors$pid_prior * priors$StreamOrderPrior * priors$PresencePrior 
      
    } else if (watershed == "Yukon") {
      gen_prior <- rep(0, length = length(pid_iso))
      gen_prior[priors$LYsites] <- as.numeric(natal_data$Lower[i])
      gen_prior[priors$MYsites] <- as.numeric(natal_data$Middle[i])
      gen_prior[priors$UYsites] <- as.numeric(natal_data$Upper[i])
      
      assign <- (1/sqrt(2*pi*error^2)) * exp(-1*(iso_o - pid_iso)^2/(2*error^2)) * 
        priors$pid_prior * priors$StreamOrderPrior * gen_prior
    }
    
    # Normalize and apply threshold
    assign_norm <- assign / sum(assign)
    assign_rescaled <- assign_norm / max(assign_norm)
    assign_rescaled[assign_rescaled < sensitivity_threshold] <- 0
    
    # Weight by CPUE
    assignment_matrix[,i] <- assign_rescaled * as.numeric(natal_data$COratio[i])
  }
  
  return(assignment_matrix)
}

#' Divide natal data into DOY quartiles (EXACT boundaries from original)
divide_doy_quartiles <- function(natal_data) {
  # EXACT quartile boundaries from original analysis
  # These match the original fixed boundaries used in the codebase
  june_11_doy <- as.numeric(format(as.Date("2020-06-11"), "%j"))  # 163
  june_21_doy <- as.numeric(format(as.Date("2020-06-21"), "%j"))  # 173  
  july_01_doy <- as.numeric(format(as.Date("2020-07-01"), "%j"))  # 183
  
  quartile_subsets <- list(
    natal_data %>% filter(DOY <= june_11_doy),                    # Q1: Start to Jun 11
    natal_data %>% filter(DOY > june_11_doy & DOY <= june_21_doy), # Q2: Jun 12-21
    natal_data %>% filter(DOY > june_21_doy & DOY <= july_01_doy), # Q3: Jun 22-Jul 1
    natal_data %>% filter(DOY > july_01_doy)                      # Q4: Jul 2-End
  )
  
  # EXACT labels from original
  labels <- c(
    sprintf("DOY Q1: Start-%d (to Jun 11)", june_11_doy),
    sprintf("DOY Q2: %d-%d (Jun 12-21)", june_11_doy + 1, june_21_doy),
    sprintf("DOY Q3: %d-%d (Jun 22-Jul 1)", june_21_doy + 1, july_01_doy),
    sprintf("DOY Q4: %d-End (Jul 2+)", july_01_doy + 1)
  )
  
  return(list(subsets = quartile_subsets, labels = labels))
}

#' Process assignments by management unit
process_mgmt_assignments <- function(edges, basin_assign) {
  if (!"mgmt_river" %in% colnames(edges)) {
    warning("No management river data found")
    return(NULL)
  }
  
  edges$basin_assign <- basin_assign
  
  # Filter to managed edges only
  managed_edges <- edges %>% 
    filter(!is.na(mgmt_river), mgmt_river != "")
  
  if (nrow(managed_edges) == 0) return(NULL)
  
  # Aggregate by management unit
  mgmt_summary <- managed_edges %>%
    st_drop_geometry() %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(basin_assign, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    ) %>%
    mutate(production_proportion = total_production / sum(total_production, na.rm = TRUE))
  
  return(mgmt_summary)
}

#' Create output directories
create_output_dirs <- function() {
  dir.create(PATHS$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$maps_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$figures_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("✓ Setup complete. All functions and parameters loaded.\n")
cat("✓ Update BASE_DIR in PATHS list for your system before running analyses.\n")
cat("✓ Remember to load 05_visualization_functions.R for plotting functions.\n")
