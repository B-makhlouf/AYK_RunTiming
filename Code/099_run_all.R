################################################################################
# 99_RUN_ALL.R - SIMPLE ANALYSIS ORCHESTRATION
################################################################################
# This script runs all analyses in sequence
# Produces identical outputs to the original complex codebase
################################################################################

cat("=== SALMON RUN TIMING ANALYSIS SUITE ===\n")

################################################################################
# SETUP CHECK
################################################################################

# Set the code directory path
CODE_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Code"

# Load all required scripts
required_scripts <- c("00_setup.R", "05_visualization_functions.R", 
                      "01_tributary_maps.R", "02_average_maps.R", 
                      "03_cumulative_distribution.R", "04_closure_protection.R")

cat("Loading required scripts...\n")
for (script in required_scripts) {
  script_path <- file.path(CODE_DIR, script)
  if (file.exists(script_path)) {
    source(script_path)
    cat(paste("✓", script, "\n"))
  } else {
    stop(paste("Missing required script:", script_path))
  }
}

# Check if setup loaded correctly
if (!exists("CONFIG") || !exists("PATHS")) {
  stop("Setup failed. Please check 00_setup.R")
}

cat("✓ All scripts loaded successfully\n\n")

################################################################################
# MAIN FUNCTION - RUN ALL ANALYSES
################################################################################

#' Run all analyses in sequence
run_all_analyses <- function(years = CONFIG$years, watersheds = CONFIG$watersheds) {
  
  cat("Starting complete analysis suite...\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n\n")
  
  start_time <- Sys.time()
  
  # Create output directories
  create_output_dirs()
  
  #--------------------------------------------------------------------------
  # 1. TRIBUTARY MAPPING ANALYSIS
  #--------------------------------------------------------------------------
  cat("1. Running tributary mapping analysis...\n")
  
  tryCatch({
    run_tributary_analysis(years, watersheds, export_csv = TRUE)
    cat("   ✓ Tributary maps and management river analysis complete\n\n")
  }, error = function(e) {
    cat("   ❌ Tributary analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # 2. AVERAGE PRODUCTION MAPS  
  #--------------------------------------------------------------------------
  cat("2. Running average production mapping analysis...\n")
  
  tryCatch({
    run_average_maps_analysis(years, watersheds)
    cat("   ✓ Average maps and boxplots complete\n\n")
  }, error = function(e) {
    cat("   ❌ Average maps analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # 3. CUMULATIVE DISTRIBUTION
  #--------------------------------------------------------------------------
  cat("3. Running cumulative distribution analysis...\n")
  
  tryCatch({
    run_cumulative_analysis(years, watersheds, CONFIG$cumulative_interval_days, export_csv = TRUE)
    cat("   ✓ Cumulative timing analysis complete\n\n")
  }, error = function(e) {
    cat("   ❌ Cumulative analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # 4. CLOSURE PROTECTION
  #--------------------------------------------------------------------------
  cat("4. Running front-end closure protection analysis...\n")
  
  tryCatch({
    run_closure_analysis(years, watersheds, CONFIG$closure_dates$start, CONFIG$closure_dates$end, export_csv = TRUE)
    cat("   ✓ Closure protection analysis complete\n\n")
  }, error = function(e) {
    cat("   ❌ Closure analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # COMPLETION SUMMARY
  #--------------------------------------------------------------------------
  end_time <- Sys.time()
  total_time <- round(difftime(end_time, start_time, units = "mins"), 1)
  
  cat("=== ANALYSIS COMPLETE ===\n")
  cat("Total time:", total_time, "minutes\n\n")
  
  cat("Check these directories for outputs:\n")
  cat("  Maps:", PATHS$maps_dir, "\n")
  cat("  Figures:", PATHS$figures_dir, "\n") 
  cat("  Data:", PATHS$output_dir, "\n\n")
  
  cat("Key outputs created:\n")
  cat("  • Tributary maps by quartile\n")
  cat("  • Management river maps and CSV data\n")
  cat("  • Average production maps and boxplots\n")
  cat("  • Cumulative timing progression plots\n")
  cat("  • Run duration analysis\n")
  cat("  • Closure protection boxplots and data\n")
}

################################################################################
# QUICK RUN FUNCTIONS
################################################################################

#' Run just one analysis
run_single_analysis <- function(analysis) {
  
  years <- CONFIG$years
  watersheds <- CONFIG$watersheds
  
  switch(analysis,
         "tributary" = run_tributary_analysis(years, watersheds, export_csv = TRUE),
         "average" = run_average_maps_analysis(years, watersheds),
         "cumulative" = run_cumulative_analysis(years, watersheds, CONFIG$cumulative_interval_days, export_csv = TRUE),
         "closure" = run_closure_analysis(years, watersheds, CONFIG$closure_dates$start, CONFIG$closure_dates$end, export_csv = TRUE),
         stop("Analysis must be one of: tributary, average, cumulative, closure")
  )
  
  cat("✓", analysis, "analysis complete\n")
}

#' Quick test run (just one year)
test_run <- function() {
  cat("Running test analysis with 2021 data only...\n")
  
  test_years <- c(2021)
  test_watersheds <- c("Kusko")
  
  cat("1. Tributary analysis...\n")
  run_tributary_analysis(test_years, test_watersheds, export_csv = TRUE)
  
  cat("2. Average maps...\n") 
  run_average_maps_analysis(test_years, test_watersheds)
  
  cat("✓ Test run complete\n")
}

################################################################################
# EXECUTION
################################################################################

cat("Available functions:\n")
cat("  run_all_analyses()        # Run complete analysis suite\n")
cat("  run_single_analysis(name) # Run one analysis: 'tributary', 'average', 'cumulative', 'closure'\n")

run_single_analysis("cumulative")
cat("  test_run()                # Quick test with 2021 data only\n\n")

cat("To run everything:\n")
cat("  run_all_analyses()\n\n")

# Uncomment to run automatically:
run_all_analyses()
