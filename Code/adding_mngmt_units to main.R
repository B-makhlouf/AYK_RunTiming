library(sf)
library(ggplot2)
library(dplyr)

# Load the shapefiles
Original <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC.shp")
W_mngmt_units <- st_read("/Users/benjaminmakhlouf/Spatial Data/Management_Units/All_Management_tribs.shp")


# Delete mgmt_river from original 
Original <- Original %>%
  select(-mgmt_river)



# Plot W_mngmt_units to visualize by mgmt_river
mgmt_units_plot <- ggplot() +
  geom_sf(data = W_mngmt_units, aes(color = mgmt_river)) +
  theme_minimal()

ggsave(
  filename = "/Users/benjaminmakhlouf/Spatial Data/Management_Units/W_mngmt_units.png",
  plot = mgmt_units_plot,
  width = 10, height = 8
)

# Ensure both layers are 2D (drop Z or M dimensions if present)
Original <- st_zm(Original, drop = TRUE, what = "ZM")
W_mngmt_units <- st_zm(W_mngmt_units, drop = TRUE, what = "ZM")

# Ensure both layers use the same CRS
if (st_crs(Original) != st_crs(W_mngmt_units)) {
  W_mngmt_units <- st_transform(W_mngmt_units, st_crs(Original))
}

# Perform spatial join - left join: keep all features from Original
Original_joined <- st_join(Original, W_mngmt_units["mgmt_river"], left = TRUE)


# make another map to make sure it worked 
mgmt_units_joined_plot <- ggplot() +
  geom_sf(data = Original_joined, aes(color = mgmt_river)) +
  theme_minimal()

# save as a png
ggsave(
  filename = "/Users/benjaminmakhlouf/Spatial Data/Management_Units/Original_joined.png",
  plot = mgmt_units_joined_plot,
  width = 10, height = 8
)


# Define output path (without extension)
output_dir <- "/Users/benjaminmakhlouf/Spatial Data"
output_name <- "KuskoUSGS_HUC_joined"
output_path <- file.path(output_dir, output_name)

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Delete old shapefile components if they exist
shp_extensions <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
for (ext in shp_extensions) {
  file_to_delete <- paste0(output_path, ext)
  if (file.exists(file_to_delete)) {
    file.remove(file_to_delete)
  }
}

# Write shapefile
st_write(Original_joined, paste0(output_path, ".shp"))
