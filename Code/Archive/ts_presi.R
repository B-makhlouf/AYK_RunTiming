### Full returns 

totalruns<- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Data/Escapements_runsize.csv")

# Plot of total Run by year 
library(ggplot2)
ggplot(totalruns, aes(x=Year, y=Total.Run)) +
  geom_line(color = "#FB8B24", size = 1.2) +
  geom_point(color = "#FB8B24", size = 2) +
  labs(title="Kuskokwim Total Run Size by Year", x="Year", y="Total Run Size") +
  theme_minimal() 
  

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





# ============================================================================
# CPUE QUARTILE TIME SERIES VISUALIZATION
# Shows seasonal timing patterns of salmon run across years
# ============================================================================

library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(RColorBrewer)

# Load the data
data_path <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
mgmt_data <- read_csv(data_path)

# ============================================================================
# DATA PREPARATION
# ============================================================================


mgmt_data<- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")




# Create summary of CPUE proportions by quartile and year
# Since cpue_prop_in_quartile should be the same for all management rivers 
# within a given year-quartile combination, we'll take the unique values
cpue_summary <- mgmt_data %>%
  select(year, quartile, cpue_prop_in_quartile) %>%
  distinct() %>%
  arrange(year, quartile)

# Check the data structure
cat("Data structure:\n")
print(head(cpue_summary))
cat("\nUnique years:", paste(sort(unique(cpue_summary$year)), collapse = ", "), "\n")
cat("Unique quartiles:", paste(sort(unique(cpue_summary$quartile)), collapse = ", "), "\n")
cat("Total observations:", nrow(cpue_summary), "\n")

# Create numeric quartile for proper ordering and continuous time variable
cpue_summary <- cpue_summary %>%
  mutate(
    quartile_num = case_when(
      quartile == "Q1" ~ 1,
      quartile == "Q2" ~ 2,
      quartile == "Q3" ~ 3,
      quartile == "Q4" ~ 4,
      TRUE ~ as.numeric(gsub("Q", "", quartile))
    ),
    # Create continuous time variable for smooth plotting
    time_continuous = year + (quartile_num - 1) * 0.25,
    # Create discrete time labels for x-axis
    time_label = paste0(year, "-", quartile),
    # Convert to percentage for easier interpretation
    cpue_percent = cpue_prop_in_quartile * 100
  ) %>%
  arrange(year, quartile_num)

# ============================================================================
# VISUALIZATION 1: SINGLE CONTINUOUS TIME SERIES LINE
# ============================================================================

# Create time labels for x-axis (every other point to avoid crowding)
time_breaks <- cpue_summary$time_continuous[seq(1, nrow(cpue_summary), by = 2)]
time_labels <- cpue_summary$time_label[seq(1, nrow(cpue_summary), by = 2)]

p1 <- ggplot(cpue_summary, aes(x = time_continuous, y = cpue_percent)) +
  geom_line(size = 1.5, color = "steelblue", alpha = 0.8) +
  geom_point(size = 3, color = "steelblue", alpha = 0.9) +
  scale_x_continuous(
    name = "Time Period",
    breaks = time_breaks,
    labels = time_labels
  ) +
  scale_y_continuous(
    name = "CPUE Proportion (%)",
    labels = scales::number_format(accuracy = 0.1, suffix = "%"),
    limits = c(0, max(cpue_summary$cpue_percent) * 1.05)
  ) +
  labs(
    title = "CPUE vs Quartile 2017-2021",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p1)
