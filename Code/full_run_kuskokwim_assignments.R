################################################################################
# KUSKOKWIM FULL RUN NATAL ORIGIN ASSIGNMENTS
################################################################################
# Single script to process all years of natal origin data and export
# scaled assignments (sum to 1) for the Kuskokwim watershed
#
# Author: Generated for Kuskokwim natal origin analysis
# Date: 2025-10-21
################################################################################

cat("=== KUSKOKWIM FULL RUN ASSIGNMENTS ===\n\n")

################################################################################
# 1. LOAD LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(glue)
})

################################################################################
# 2. SET PATHS AND PARAMETERS
################################################################################

# Base directory - UPDATE THIS FOR YOUR SYSTEM
BASE_DIR <- "/Users/benjaminmakhlouf"

# Data paths
KUSKO_EDGES_PATH <- file.path(BASE_DIR, "Spatial Data/KuskoUSGS_HUC_joined.shp")
KUSKO_BASIN_PATH <- file.path(BASE_DIR, "Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp")
NATAL_DATA_DIR <- file.path(BASE_DIR,
                             "Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/Data",
                             "Natal Origin Analysis Data/03_Natal Origins Genetics CPUE")

# Output directory
OUTPUT_DIR <- "Analysis_Results/Full_Run_Assignments"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
YEARS <- c(2017, 2018, 2019, 2020, 2021, 2022)
SENSITIVITY_THRESHOLD <- 0.7  # 70% threshold for assignments
MIN_ERROR <- 0.0006           # Minimum isotope error
MIN_STREAM_ORDER <- 3         # Minimum stream order to include

# Error calculation constants
WITHIN_SITE_ERROR <- 0.0003133684 / 1.96
ANALYTICAL_ERROR <- 0.00011 / 2

################################################################################
# 3. LOAD SPATIAL DATA (ONCE)
################################################################################

cat("Loading spatial data...\n")

# Load edges and basin
edges <- st_read(KUSKO_EDGES_PATH, quiet = TRUE)
basin <- st_read(KUSKO_BASIN_PATH, quiet = TRUE)

# Transform and filter by stream order
edges <- st_transform(edges, st_crs(basin))
edges <- edges[edges$Str_Order >= MIN_STREAM_ORDER, ]

# Consolidate Johnson into Lower Kusko (management unit consolidation)
edges$mgmt_river[edges$mgmt_river == "Johnson"] <- "Lower Kusko"

# Extract isoscape predictions for assignments
pid_iso <- edges$iso_pred
pid_isose <- edges$isose_pred

cat(glue("  Loaded {nrow(edges)} stream edges\n"))
cat(glue("  Unique management units: {length(unique(edges$mgmt_river[!is.na(edges$mgmt_river)]))}\n\n"))

################################################################################
# 4. SETUP PRIORS (ONCE)
################################################################################

cat("Setting up spatial priors...\n")

# Stream order prior: Include only streams >= min stream order
StreamOrderPrior <- ifelse(edges$Str_Order >= MIN_STREAM_ORDER, 1, 0)

# Isotope prior: Use UniPh2oNoE for Kuskokwim
pid_prior <- edges$UniPh2oNoE

# Presence prior: Exclude high-order streams without spawning capacity
PresencePrior <- ifelse((edges$Str_Order %in% c(6, 7, 8)) & edges$SPAWNING_C == 0, 0, 1)

# Calculate isotope error for each edge
pid_isose_mod <- ifelse(pid_isose < MIN_ERROR, MIN_ERROR, pid_isose)
error <- sqrt(pid_isose_mod^2 + WITHIN_SITE_ERROR^2 + ANALYTICAL_ERROR^2)

cat("  Priors configured\n\n")

################################################################################
# 5. STORAGE FOR RESULTS
################################################################################

# Initialize lists to store results
all_edge_assignments <- list()
all_mgmt_assignments <- list()

################################################################################
# 6. MAIN PROCESSING LOOP - ALL YEARS
################################################################################

cat("Processing natal origin assignments for all years...\n\n")

for (year in YEARS) {

  cat(glue("--- YEAR {year} ---\n"))

  # ============================================================================
  # 6.1 LOAD NATAL DATA FOR THIS YEAR
  # ============================================================================

  natal_file <- file.path(NATAL_DATA_DIR, paste0(year, "_Kusko_Natal_Origins_Genetics_CPUE.csv"))

  if (!file.exists(natal_file)) {
    cat(glue("  WARNING: Natal data not found for {year}, skipping\n\n"))
    next
  }

  natal_data <- read_csv(natal_file, show_col_types = FALSE)

  # Clean data: Remove fish with missing isotope or CPUE data
  natal_data <- natal_data %>%
    filter(!is.na(natal_iso), !is.na(dailyCPUEprop), !is.na(COratio))

  n_fish <- nrow(natal_data)
  cat(glue("  Loaded {n_fish} fish observations\n"))

  if (n_fish == 0) {
    cat("  No valid fish data, skipping\n\n")
    next
  }

  # ============================================================================
  # 6.2 PERFORM BAYESIAN ASSIGNMENTS FOR ALL FISH
  # ============================================================================

  cat("  Performing Bayesian assignments...\n")

  # Initialize assignment matrix: rows = stream edges, columns = individual fish
  assignment_matrix <- matrix(NA, nrow = length(pid_iso), ncol = n_fish)

  # Loop through each fish and calculate assignment probabilities
  for (i in 1:n_fish) {

    # Get this fish's isotope value
    iso_o <- as.numeric(natal_data$natal_iso[i])

    # BAYESIAN ASSIGNMENT CALCULATION
    # Combines: isotope likelihood + spatial priors
    assign <- (1/sqrt(2*pi*error^2)) *
              exp(-1*(iso_o - pid_iso)^2/(2*error^2)) *
              pid_prior *
              StreamOrderPrior *
              PresencePrior

    # Normalize to sum to 1
    assign_norm <- assign / sum(assign)

    # Rescale so maximum = 1
    assign_rescaled <- assign_norm / max(assign_norm)

    # Apply sensitivity threshold (set values < 0.7 to 0)
    assign_rescaled[assign_rescaled < SENSITIVITY_THRESHOLD] <- 0

    # Weight by catch-origin ratio (CPUE weighting)
    assignment_matrix[,i] <- assign_rescaled * as.numeric(natal_data$COratio[i])
  }

  # ============================================================================
  # 6.3 AGGREGATE ASSIGNMENTS ACROSS ALL FISH
  # ============================================================================

  cat("  Aggregating across fish...\n")

  # Sum assignments across all fish for each stream edge
  basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)

  # Scale to proportions (sum to 1) - THIS IS THE KEY OUTPUT
  basin_assign_rescale <- basin_assign_sum / sum(basin_assign_sum, na.rm = TRUE)

  # Also calculate normalized version (max = 1) for visualization
  basin_assign_norm <- basin_assign_rescale / max(basin_assign_rescale, na.rm = TRUE)

  # ============================================================================
  # 6.4 STORE EDGE-LEVEL RESULTS
  # ============================================================================

  # Create data frame with edge-level assignments
  edge_results <- data.frame(
    year = year,
    edge_id = 1:length(basin_assign_sum),
    Str_Order = edges$Str_Order,
    mgmt_river = edges$mgmt_river,
    iso_pred = edges$iso_pred,
    assignment_sum = basin_assign_sum,
    assignment_proportion = basin_assign_rescale,
    assignment_normalized = basin_assign_norm
  )

  all_edge_assignments[[as.character(year)]] <- edge_results

  # ============================================================================
  # 6.5 AGGREGATE TO MANAGEMENT UNITS
  # ============================================================================

  cat("  Aggregating to management units...\n")

  # Add assignments to edges spatial data
  edges_temp <- edges
  edges_temp$basin_assign <- basin_assign_rescale

  # Filter to managed edges only
  managed_edges <- edges_temp %>%
    filter(!is.na(mgmt_river), mgmt_river != "")

  # Aggregate by management unit
  mgmt_summary <- managed_edges %>%
    st_drop_geometry() %>%
    group_by(mgmt_river) %>%
    summarise(
      total_production = sum(basin_assign, na.rm = TRUE),
      edge_count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      production_proportion = total_production / sum(total_production, na.rm = TRUE),
      year = year
    ) %>%
    select(year, mgmt_river, production_proportion, total_production, edge_count)

  all_mgmt_assignments[[as.character(year)]] <- mgmt_summary

  cat(glue("  Processed {nrow(mgmt_summary)} management units\n"))
  cat(glue("  Total production proportion sum: {round(sum(mgmt_summary$production_proportion), 4)}\n\n"))

}

################################################################################
# 7. COMBINE AND EXPORT RESULTS
################################################################################

cat("Exporting results...\n")

# Combine all years for edge-level assignments
edge_assignments_all <- bind_rows(all_edge_assignments)

# Combine all years for management unit assignments
mgmt_assignments_all <- bind_rows(all_mgmt_assignments)

# Export edge-level assignments
edge_output_file <- file.path(OUTPUT_DIR, "kuskokwim_edge_assignments_full_run.csv")
write_csv(edge_assignments_all, edge_output_file)
cat(glue("  Saved edge-level assignments: {edge_output_file}\n"))

# Export management unit assignments
mgmt_output_file <- file.path(OUTPUT_DIR, "kuskokwim_management_unit_assignments_full_run.csv")
write_csv(mgmt_assignments_all, mgmt_output_file)
cat(glue("  Saved management unit assignments: {mgmt_output_file}\n"))

################################################################################
# 8. SUMMARY STATISTICS
################################################################################

cat("\n=== SUMMARY ===\n")
cat(glue("Years processed: {paste(unique(edge_assignments_all$year), collapse=', ')}\n"))
cat(glue("Total observations: {nrow(edge_assignments_all)}\n"))
cat(glue("Stream edges per year: {nrow(edge_assignments_all) / length(YEARS)}\n"))
cat(glue("Management units: {length(unique(mgmt_assignments_all$mgmt_river))}\n"))
cat("\nManagement units (upstream to downstream):\n")
for (unit in unique(mgmt_assignments_all$mgmt_river)) {
  cat(glue("  - {unit}\n"))
}

cat("\n=== COMPLETE ===\n")
cat("Scaled assignments (sum to 1) exported successfully.\n")
