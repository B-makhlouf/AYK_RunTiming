library(sdmTMB)
library(dplyr)
library(sf)
library(tidyr)
library(ggplot2)
library(viridisLite)
library(viridis)

## Read in the Kusko edges shapefile 
kusko_edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/Cleaned AYK Shapefiles/Kusko_cleaned.shp")
## Read in the Basin 
kusko_basin <- st_read("/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp")

# Load in all production data 
ALL_prod_Data <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/all_doy_tributary_values.csv")

# Function to extract points on lines for a given set of edges
# Updated to ensure rid is preserved
extract_points_on_lines <- function(edges) {
  edges_with_points <- edges %>%
    mutate(point_on_line = st_point_on_surface(geometry))
  
  point_coords <- st_coordinates(edges_with_points$point_on_line)
  colnames(point_coords) <- c("X", "Y")
  
  edges_with_coords <- edges_with_points %>%
    st_drop_geometry() %>%
    bind_cols(as.data.frame(point_coords))
  
  return(edges_with_coords)
}

# Create a list to store each year/quartile combination
tidy_data_list <- list()

# Define the years and quartiles to process
years_to_process <- c(2017, 2018, 2019, 2020, 2021)
quartiles_to_process <- c("Q1", "Q2", "Q3", "Q4")
numeric_quartiles <- c(1, 2, 3, 4)  # Corresponding numeric values

# Create a mapping between string quartiles and numeric quartiles
quartile_mapping <- setNames(numeric_quartiles, quartiles_to_process)

# Process each year and quartile combination
for (i in 1:length(quartiles_to_process)) {
  quartile_val <- quartiles_to_process[i]
  numeric_quartile <- numeric_quartiles[i]
  
  for (year_val in years_to_process) {
    # Check if year exists in the data
    if (!year_val %in% ALL_prod_Data$year) {
      message(paste("Year", year_val, "not found in data, skipping."))
      next
    }
    
    # Filter data for this year and quartile
    current_prod_data <- ALL_prod_Data %>% 
      filter(year == year_val) %>% 
      filter(quartile == quartile_val) %>% 
      filter(analysis_type == "Cumulative_DOY")
    
    # Skip if no data for this combination
    if (nrow(current_prod_data) == 0) {
      message(paste("No data for Year", year_val, "Quartile", quartile_val, ", skipping."))
      next
    }
    
    # Filter edges to match the stream order in the current data
    min_strOrd <- min(current_prod_data$stream_order, na.rm = TRUE)
    current_edges <- kusko_edges %>% filter(Str_Order >= min_strOrd)
    
    # Check if lengths match - use appropriate joining if they don't
    if (length(current_edges$Str_Order) != length(current_prod_data$stream_order)) {
      message(paste("Length mismatch for Year", year_val, "Quartile", quartile_val))
      
      # Try to match by segment ID if available
      if ("segment_id" %in% colnames(current_prod_data) && "segment_id" %in% colnames(current_edges)) {
        message("Attempting to join by segment_id...")
      } else {
        warning(paste("Cannot match data for Year", year_val, "Quartile", quartile_val))
        next
      }
    }
    
    # Extract points on the lines
    edges_with_coords <- extract_points_on_lines(current_edges)
    
    # Create data frame for this year/quartile combination
    # Now including rid as a unique identifier
    current_tidy_data <- data.frame(
      X = edges_with_coords$X,
      Y = edges_with_coords$Y,
      X_km = edges_with_coords$X / 1000,
      Y_km = edges_with_coords$Y / 1000,
      rid = edges_with_coords$rid,  # Add the rid here
      Year = year_val,
      Quartile = numeric_quartile,  # Using the numeric value
      prod_value = current_prod_data$basin_assign_norm,
      segment_id = ifelse("segment_id" %in% colnames(current_prod_data), 
                          current_prod_data$segment_id, NA)
    )
    
    # Add to the list
    tidy_data_list[[paste(year_val, numeric_quartile, sep = "_")]] <- current_tidy_data
  }
}

# Combine all data frames into a single tidy data frame
tidy_prod_data <- bind_rows(tidy_data_list)

# Make sure Quartile is numeric
tidy_prod_data$Quartile <- as.numeric(tidy_prod_data$Quartile)

# Print a count of unique rid values to verify
print(paste("Number of unique rid values in tidy data:", length(unique(tidy_prod_data$rid))))

################ 
################ Data is now in the correct form 

#### Make the mesh 
mesh <- make_mesh(tidy_prod_data, c("X_km", "Y_km"), cutoff = 10) 
plot(mesh)

# Fit a Tweedie spatial random field GLMM
fit <- sdmTMB(
  prod_value ~ 1,
  data = tidy_prod_data, mesh = mesh,
  family = tweedie(link = "log"),
  time = "Year",
  spatiotemporal = "iid"
)

# Create a prediction grid that includes rid
# First, get unique combinations of X, Y, and rid from the original data
unique_locations <- tidy_prod_data %>%
  select(X, Y, X_km, Y_km, rid) %>%
  distinct()

# Create replicated data frame with all years for each location
nd <- unique_locations %>%
  tidyr::crossing(Year = unique(tidy_prod_data$Year))

# Make predictions
predictions <- predict(fit, newdata = nd)

# Summarize predictions by rid to maintain all sites
pred <- predictions %>%
  group_by(rid) %>%
  summarise(
    X = first(X),
    Y = first(Y),
    sd_est = sd(est),
    cv_est = sd(est) / mean(est),
    sd_rf = sd(est_rf)
  ) %>%
  ungroup()

# Verify counts
print(paste("Number of unique locations:", nrow(unique_locations)))
print(paste("Number of sites after summarizing:", nrow(pred)))

# Create the plot with basin boundary and river network
basin_df <- st_as_sf(kusko_basin)


# Add the pred values back to kusko_edges, filtered to the min_stream_order 

kusko_edges <- kusko_edges %>%
  filter(Str_Order >= min_strOrd) %>%
  left_join(pred, by = "rid")


ggplot() +
  # Add the basin outline
  geom_sf(data = basin_df, fill = NA, color = "black", size = 0.5) +
  # Add the river network colored by standard deviation
  geom_sf(data = kusko_edges, aes(color = sd_est), size = 0.7) +
  # Use viridis color scale
  scale_color_viridis(name = "Standard\nDeviation", 
                      option = "plasma", 
                      direction = -1) +  # Reverse direction for better contrast
  # Add a title and proper labels
  labs(title = "Variability in Fish Run Timing Across the Kuskokwim River Network",
       subtitle = paste("Based on data from", min(years_to_process), "to", max(years_to_process)),
       x = "Longitude", y = "Latitude") +
  # Use a cleaner theme for maps
  theme_minimal() +
  # Additional theme adjustments for maps
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.height = unit(1, "cm")
  )





# Alternative plot using cv_est (coefficient of variation) instead of sd_est
ggplot() +
  geom_sf(data = basin_df, fill = NA, color = "black", size = 0.5) +
  geom_sf(data = kusko_edges, color = "blue", size = 0.3, alpha = 0.6) +
  geom_point(data = pred, aes(X, Y, color = cv_est)) +
  scale_color_viridis(name = "CV of\nEstimates") +
  labs(title = "Coefficient of Variation in Fish Run Timing",
       subtitle = "Higher values indicate more variable timing across years",
       x = "Longitude", y = "Latitude") +
  theme_minimal()