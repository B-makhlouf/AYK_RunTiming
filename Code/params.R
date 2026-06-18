################################################################################
# SHARED ANALYSIS PARAMETERS
#
# Single source of truth for all production map and contour map scripts.
# Edit values here — changes will propagate automatically when each script runs.
#
# Sourced by:
#   00_All_FullRun_ProductionMaps.R
#   00_ProductionMaps50pct.R
#   00_ProductionMapsQuartiles.R
#   Analysis/01_CONTOURMAPS_FullYrBothBasins.R
#   Analysis/01_CONTOURMAPS_FullYr_Kusko.R
#   Analysis/01_CONTOURMAPS_HalfYrBothBasins.R
################################################################################

# ==============================================================================
# KUSKOKWIM PARAMETERS
# ==============================================================================
KUSKO_PARAMS <- list(
  min_stream_order      = 3,       # Minimum Strahler stream order included
  min_error             = 0.0006,  # Lower bound clamp on pid_isose error
  max_error             = 0.00089, # Upper bound clamp (used in Quartiles analysis)
  sensitivity_threshold = 0.0,     # Assignment rescaled values below this -> 0
  channel_slope_cutoff  = 2.0      # NewHabitatPrior: Channel_sl > this -> excluded
)

# ==============================================================================
# YUKON PARAMETERS
# ==============================================================================
YUKON_PARAMS <- list(
  min_stream_order      = 3,       # Minimum Strahler stream order included
  min_error             = 0.0035,  # Lower bound clamp on pid_isose error
  sensitivity_threshold = 0.7,     # Assignment rescaled values below this -> 0
  channel_slope_cutoff  = 2.5,     # NewHabitatPrior: Channel_sl > this -> excluded
  porcupine_penalty     = 0.0      # Post-hoc downweight multiplier for Porc_off == 0 segments
)

# ==============================================================================
# SHARED ANALYSIS PARAMETERS
# ==============================================================================
PRODUCTION_THRESHOLD <- 0.7   # Minimum normalised production value to include a segment
TEMP_INTERVAL_DAYS   <- 3     # Temporal sampling interval for temperature extraction (days)
