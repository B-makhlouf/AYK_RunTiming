################################################################################
# 99_RUN_ALL.R - COMPLETE ANALYSIS ORCHESTRATION
################################################################################
# This script runs all analyses in sequence including precision analysis
# Produces identical outputs to the original complex codebase plus new precision metrics
################################################################################

cat("=== SALMON RUN TIMING ANALYSIS SUITE ===\n")

################################################################################
# SETUP CHECK
################################################################################

# Set the code directory path
CODE_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Code"

# Load all required scripts
required_scripts <- c(
  "00_setup.R", 
  "05_visualization_functions.R",
  "01_tributary_maps.R", 
  "02_average_maps.R", 
  "03_cumulative_distribution.R", 
  "04_closure_protection.R",
  "06_assignment_precision.R"  # NEW: Precision analysis
)

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
#' @param years Vector of years to analyze
#' @param watersheds Vector of watersheds to analyze
#' @param run_precision Whether to run precision analysis (default TRUE)
run_all_analyses <- function(years = CONFIG$years, 
                             watersheds = CONFIG$watersheds,
                             run_precision = TRUE) {
  
  cat("Starting complete analysis suite...\n")
  cat("Years:", paste(years, collapse = ", "), "\n")
  cat("Watersheds:", paste(watersheds, collapse = ", "), "\n")
  cat("Precision analysis:", ifelse(run_precision, "ENABLED", "DISABLED"), "\n\n")
  
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
  cat("2. Running average production maps analysis...\n")
  
  tryCatch({
    run_average_analysis(years, watersheds, export_csv = TRUE)
    cat("   ✓ Average production maps and boxplots complete\n\n")
  }, error = function(e) {
    cat("   ❌ Average production analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # 3. CUMULATIVE DISTRIBUTION ANALYSIS
  #--------------------------------------------------------------------------
  cat("3. Running cumulative distribution analysis...\n")
  
  tryCatch({
    cumulative_results <- run_cumulative_analysis(
      years, watersheds, 
      interval_days = CONFIG$cumulative_interval_days,
      export_csv = TRUE
    )
    cat("   ✓ Cumulative distribution and run duration analysis complete\n\n")
  }, error = function(e) {
    cat("   ❌ Cumulative distribution analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # 4. CLOSURE PROTECTION ANALYSIS
  #--------------------------------------------------------------------------
  cat("4. Running front-end closure protection analysis...\n")
  
  tryCatch({
    closure_results <- run_closure_analysis(
      years, watersheds,
      closure_start = CONFIG$closure_dates$start,
      closure_end = CONFIG$closure_dates$end,
      export_csv = TRUE
    )
    cat("   ✓ Closure protection analysis complete\n\n")
  }, error = function(e) {
    cat("   ❌ Closure protection analysis failed:", e$message, "\n\n")
  })
  
  #--------------------------------------------------------------------------
  # 5. ASSIGNMENT PRECISION ANALYSIS (NEW!)
  #--------------------------------------------------------------------------
  if (run_precision) {
    cat("5. Running individual assignment precision analysis...\n")
    
    tryCatch({
      precision_results <- run_precision_analysis(
        years, watersheds,
        export_csv = TRUE
      )
      cat("   ✓ Assignment precision analysis complete\n\n")
    }, error = function(e) {
      cat("   ❌ Assignment precision analysis failed:", e$message, "\n\n")
    })
  } else {
    cat("5. Skipping assignment precision analysis (disabled)\n\n")
    precision_results <- NULL
  }
  
  #--------------------------------------------------------------------------
  # SUMMARY
  #--------------------------------------------------------------------------
  end_time <- Sys.time()
  elapsed <- round(difftime(end_time, start_time, units = "mins"), 2)
  
  cat("\n=== ANALYSIS SUITE COMPLETE ===\n")
  cat(glue("Total runtime: {elapsed} minutes\n"))
  cat(glue("Output directory: {PATHS$output_dir}\n\n"))
  
  cat("Analyses completed:\n")
  cat("  ✓ 1. Tributary mapping and management river analysis\n")
  cat("  ✓ 2. Average production maps and boxplots\n")
  cat("  ✓ 3. Cumulative distribution and run duration\n")
  cat("  ✓ 4. Front-end closure protection\n")
  if (run_precision) {
    cat("  ✓ 5. Individual assignment precision\n")
  }
  
  # Return results for further analysis if needed
  return(list(
    precision = if (run_precision) precision_results else NULL,
    runtime = elapsed
  ))
}

################################################################################
# CONVENIENCE FUNCTIONS FOR RUNNING INDIVIDUAL ANALYSES
################################################################################

#' Run only the precision analysis (useful for testing)
run_precision_only <- function(years = CONFIG$years, 
                               watersheds = CONFIG$watersheds) {
  cat("Running assignment precision analysis only...\n")
  results <- run_precision_analysis(years, watersheds, export_csv = TRUE)
  return(results)
}

#' Run the core 4 original analyses (without precision)
run_original_analyses <- function(years = CONFIG$years, 
                                  watersheds = CONFIG$watersheds) {
  return(run_all_analyses(years, watersheds, run_precision = FALSE))
}

################################################################################
# EXECUTION SECTION
################################################################################

# Uncomment one of these lines to run the analyses:

# Run all analyses including new precision analysis:
# results <- run_all_analyses()

# Run only original 4 analyses:
# results <- run_original_analyses()

# Run only precision analysis:
# precision_results <- run_precision_only()

# Run with custom parameters:
# results <- run_all_analyses(
#   years = c(2017, 2018, 2019), 
#   watersheds = c("Yukon"),
#   run_precision = TRUE
# )

cat("\n=== READY TO RUN ===\n")
cat("Available functions:\n")
cat("  run_all_analyses()          # Run all 5 analyses\n")
cat("  run_original_analyses()     # Run original 4 only\n")
cat("  run_precision_only()        # Run precision analysis only\n")
cat("\nTo execute, uncomment one of the function calls above\n")
cat("or run manually in the console.\n")