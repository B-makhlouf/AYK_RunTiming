################################################################################
# SIMPLE SHAPEFILE EXPORT WITH AVERAGE PRODUCTION VALUES
################################################################################

library(sf)
library(dplyr)
library(readr)
library(tidyr)

# Load data
edges <- st_read("/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp", quiet = TRUE)
avg_data <- read_csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Average_Management_Production/2017-2022/average_production_by_quartile_2017-2022.csv", show_col_types = FALSE)

# Reshape data to wide format
avg_wide <- avg_data %>%
  select(mgmt_river, quartile_clean, avg_within_quartile_prop) %>%
  pivot_wider(names_from = quartile_clean, values_from = avg_within_quartile_prop, names_prefix = "avg_") %>%
  mutate_at(vars(starts_with("avg_")), ~round(coalesce(., 0), 4))

# Join and save
result <- edges %>% left_join(avg_wide, by = "mgmt_river")
st_write(result, "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/kusko_with_averages.shp", delete_dsn = TRUE)

cat("✓ Done! Saved to: kusko_with_averages.shp\n")