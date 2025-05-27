library(sf)

# Load the original edges shapefile
Original <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC.shp")
W_mngmt_units <- st_read("/Users/benjaminmakhlouf/Spatial Data/Management_Units/kusko_edges_with_management.shp")

# Add mgmt_river to Original from W_mngmt_units
Original$mgmt_river <- W_mngmt_units$mgmt_river

# Convert from 3D to 2D geometry
Original <- st_zm(Original)

# Check if directory exists and create if needed
output_dir <- "/Users/benjaminmakhlouf/Spatial Data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Use file.path for better path handling
output_path <- file.path(output_dir, "KuskoUSGS_HUC.shp")

# Overwrite the original file
st_write(Original, output_path, delete_dsn = TRUE)
