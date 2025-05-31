### Full returns 

totalruns<- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Escapements_runsize.csv")

# Plot of total Run by year 
library(ggplot2)
ggplot(totalruns, aes(x=Year, y=Total.Run)) +
  geom_line(color = "#FB8B24", size = 1.2) +
  geom_point(color = "#FB8B24", size = 2) +
  labs(title="Kuskokwim Total Run Size by Year", x="Year", y="Total Run Size") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#2E3440", color = NA),
    panel.background = element_rect(fill = "#2E3440", color = NA),
    panel.grid = element_line(color = "#3B4252"),
    text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    plot.title = element_text(color = "white")
  )


# annual_total_run_maps.R
# Create annual basin maps showing tributaries colored by total run size
# Simple visualization of how total run size changes across years

library(sf)
library(dplyr)
library(here)
library(RColorBrewer)

# Source required utilities (just for loading spatial data)
source(here("code/utils/spatial_utils.R"))

#' Create annual basin maps colored by total run size
#' @param years Vector of years to process
#' @param watershed Watershed name ("Kusko" or "Yukon")
#' @return Invisibly returns the run size data
create_annual_total_run_maps <- function(years = c(2017, 2018, 2019, 2020, 2021), 
                                         watershed = "Kusko") {
  
  message("=== Creating Annual Total Run Maps ===")
  
  # Load total run data
  totalruns <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Escapements_runsize.csv")
  
  # Filter for the watershed and years
  run_data <- totalruns %>% 
    filter(River == watershed, Year %in% years)
  
  if (nrow(run_data) == 0) {
    stop(paste("No run size data found for", watershed, "in years:", paste(years, collapse = ", ")))
  }
  
  print("Run size data:")
  print(run_data)
  
  # Get watershed parameters (for stream order filtering)
  if (watershed == "Kusko") {
    min_stream_order <- 3
  } else if (watershed == "Yukon") {
    min_stream_order <- 5
  } else {
    stop("Watershed must be 'Kusko' or 'Yukon'")
  }
  
  # Create output directory
  output_dir <- here("Basin Maps/Annual_Total_Run")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load spatial data (same for all years)
  spatial_data <- load_spatial_data(watershed, 8, min_stream_order)
  edges <- spatial_data$edges
  basin <- spatial_data$basin
  
  # Set fixed scale for easier readability
  global_min <- 100000
  global_max <- 250000
  
  message(paste("Using fixed scale:", format(global_min, big.mark = ","), "to", format(global_max, big.mark = ",")))
  
  # Color palette (same as tributary maps)
  pallete <- brewer.pal(9, "YlOrRd")
  pallete_expanded <- colorRampPalette(pallete)(10)
  
  # ============================================================================
  # Create maps for each year
  # ============================================================================
  
  for (year in years) {
    message(paste("Creating map for", year))
    
    # Get total run for this year
    year_run <- run_data$Total.Run[run_data$Year == year]
    
    if (length(year_run) == 0) {
      message(paste("No run data for", year, "- skipping"))
      next
    }
    
    # Normalize run size to 0-1 scale using global min/max
    run_normalized <- (year_run - global_min) / (global_max - global_min)
    
    # Color ALL tributaries the same color based on that year's total run
    # Scale across years: higher total run = deeper red, lower = lighter yellow
    
    # Map normalized value (0-1) to color palette index (1-10)
    color_index <- pmax(1, pmin(10, round(run_normalized * 9) + 1))
    
    # Get the color for this year's total run
    tributary_color <- pallete_expanded[color_index]
    
    message(paste("Year", year, "- Total Run:", format(year_run, big.mark = ","), 
                  "- Normalized:", round(run_normalized, 3), 
                  "- Color Index:", color_index))
    
    # Create color coding - all tributaries same color within year
    colcode <- rep(tributary_color, nrow(edges))
    
    # Set excluded areas to gray (low stream order)
    colcode[edges$Str_Order < min_stream_order] <- 'gray60'
    
    # Set line widths by stream order (same as original)
    stream_order_lwd <- edges$Str_Order
    linewidths <- rep(1, length(stream_order_lwd))
    
    if (watershed == "Yukon") {
      linewidths <- ifelse(stream_order_lwd == 9, 3.7, linewidths)
      linewidths <- ifelse(stream_order_lwd == 8, 2.5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 7, 1.7, linewidths)
      linewidths <- ifelse(stream_order_lwd == 6, 1.5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 5, 1, linewidths)
      linewidths <- ifelse(stream_order_lwd == 4, 1, linewidths)
      linewidths <- ifelse(stream_order_lwd == 3, 1, linewidths)
    } else {
      linewidths <- ifelse(stream_order_lwd == 9, 5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 8, 4, linewidths)
      linewidths <- ifelse(stream_order_lwd == 7, 3, linewidths)
      linewidths <- ifelse(stream_order_lwd == 6, 2, linewidths)
      linewidths <- ifelse(stream_order_lwd == 5, 1.8, linewidths)
      linewidths <- ifelse(stream_order_lwd == 4, 1.5, linewidths)
      linewidths <- ifelse(stream_order_lwd == 3, 1, linewidths)
    }
    
    # Create map
    map_filename <- paste0(watershed, "_", year, "_Total_Run_", 
                           format(year_run, scientific = FALSE), ".png")
    map_path <- file.path(output_dir, map_filename)
    
    png(file = map_path, width = 12, height = 10, units = "in", res = 300, bg = "transparent")
    
    # Simple title - just the run size
    plot_title <- format(year_run, big.mark = ",")
    
    # Set plot margins - reduce top margin to bring title closer
    par(mar = c(2, 1, 2, 1), bg = "transparent")
    
    # Plot basin and tributaries - title closer to figure, transparent background
    plot(st_geometry(basin), col = 'gray60', border = 'gray60', 
         main = plot_title, bg = "transparent", cex.main = 3)
    plot(st_geometry(edges), col = colcode, pch = 16, axes = FALSE, 
         add = TRUE, lwd = linewidths)
    
    dev.off()
    
    # Reset par
    par(mar = c(5, 4, 4, 2) + 0.1, bg = "white")
    
    message(paste("Created:", map_filename))
  }
  
  message("=== Annual Total Run Maps Complete ===")
  message(paste("Maps saved to:", output_dir))
  
  # Return run data
  invisible(run_data)
}

# ============================================================================
# EXECUTE: Create the maps
# ============================================================================

# Create annual total run maps for Kuskokwim 2017-2021
run_data <- create_annual_total_run_maps(
  years = c(2017, 2018, 2019, 2020, 2021),
  watershed = "Kusko"
)

message("\n=== Summary ===")
message("Created basin maps where:")
message("- All tributaries within each year are colored identically")
message("- Color intensity reflects that year's total run size RELATIVE to other years") 
message("- Higher total run years (like 2019) = deeper red/orange")
message("- Lower total run years = lighter yellow/pale colors")
message("- Color scale consistent across all years for direct comparison")

# Print the actual color scaling for reference
message("\n=== Color Scaling by Year ===")
for (year in years) {
  year_run <- run_data$Total.Run[run_data$Year == year]
  if (length(year_run) > 0) {
    run_normalized <- (year_run - global_min) / (global_max - global_min)
    color_index <- pmax(1, pmin(10, round(run_normalized * 9) + 1))
    message(paste("Year", year, ":", format(year_run, big.mark = ","), 
                  "fish - Color intensity:", round(run_normalized * 100, 1), "%"))
  }
}