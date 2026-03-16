################################################################################
# MANAGEMENT UNIT AVERAGE ISOTOPE RATIO ANALYSIS
################################################################################
# PURPOSE: Calculate and visualize average iso_pred values by management unit
# WITH PRIORS APPLIED (only streams used in assignment)
################################################################################

library(sf)
library(dplyr)
library(ggplot2)
library(viridis)

################################################################################
# CONFIGURATION
################################################################################

# File paths - UPDATE FOR YOUR SYSTEM
EDGES_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
OUTPUT_DIR <- "Analysis_Results"

# Management unit ordering (upstream to downstream)
MGMT_ORDER <- c(
  "N. Fork Kusko", 
  "E. Fork Kuskokwim",
  "S. Fork Kusko",
  "Upper Kusko Main", 
  "Big River", 
  "Takotna and Nixon Fork",
  "Tatlawiksuk", 
  "Swift", 
  "Stony", 
  "Holitna",
  "Hoholitna",
  "Middle Kusko Main", 
  "George", 
  "Oskakawlik", 
  "Holokuk",
  "Aniak", 
  "Tuluksak", 
  "Kisaralik", 
  "Kwethluk", 
  "Lower Kusko"
)

################################################################################
# LOAD AND PROCESS DATA
################################################################################

cat("Loading spatial data...\n")
edges <- st_read(EDGES_PATH, quiet = TRUE)

# Apply priors to filter edges (mimicking assignment logic)
cat("Applying priors to filter edges...\n")

# Watershed parameters for Kusko
min_stream_order <- 3

# Setup priors (Kuskokwim-specific)
StreamOrderPrior <- ifelse(edges$Str_Order >= min_stream_order, 1, 0)
pid_prior <- edges$UniPh2oNoE
PresencePrior <- ifelse((edges$Str_Order %in% c(6, 7, 8)) & edges$SPAWNING_C == 0, 0, 1)

# Filter to edges that pass ALL prior checks
edges_with_priors <- edges %>%
  mutate(
    StreamOrderPrior = StreamOrderPrior,
    pid_prior = pid_prior,
    PresencePrior = PresencePrior
  ) %>%
  filter(
    StreamOrderPrior == 1,
    pid_prior > 0,
    PresencePrior == 1
  )

cat(paste("After applying priors:", nrow(edges_with_priors), "of", nrow(edges), "edges remain\n"))

# Filter to edges with management unit assignments
managed_edges <- edges_with_priors %>%
  filter(!is.na(mgmt_river), mgmt_river != "", !is.na(iso_pred))

cat(paste("Found", nrow(managed_edges), "stream segments with management units after priors\n"))

# Calculate average iso_pred by management unit
mgmt_iso_summary <- managed_edges %>%
  st_drop_geometry() %>%
  group_by(mgmt_river) %>%
  summarise(
    mean_iso_pred = mean(iso_pred, na.rm = TRUE),
    sd_iso_pred = sd(iso_pred, na.rm = TRUE),
    min_iso_pred = min(iso_pred, na.rm = TRUE),
    max_iso_pred = max(iso_pred, na.rm = TRUE),
    n_segments = n(),
    .groups = "drop"
  )

# Order by watershed position (upstream to downstream)
mgmt_iso_summary$mgmt_river <- factor(
  mgmt_iso_summary$mgmt_river,
  levels = MGMT_ORDER[MGMT_ORDER %in% mgmt_iso_summary$mgmt_river]
)

# Sort by watershed position
mgmt_iso_summary <- mgmt_iso_summary %>%
  arrange(mgmt_river)

################################################################################
# CREATE VISUALIZATION
################################################################################

cat("Creating visualization...\n")

# Prepare data for boxplot (need individual segment values)
boxplot_data <- managed_edges %>%
  st_drop_geometry() %>%
  filter(!is.na(mgmt_river), mgmt_river != "", !is.na(iso_pred))

# Apply management unit ordering to boxplot data (upstream at top)
boxplot_data$mgmt_river <- factor(
  boxplot_data$mgmt_river,
  levels = rev(MGMT_ORDER[MGMT_ORDER %in% boxplot_data$mgmt_river])
)

# Create color gradient matching other figures (red upstream to blue downstream)
n_units <- length(levels(boxplot_data$mgmt_river))
watershed_colors <- colorRampPalette(c(
  "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
  "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
))(n_units)
# Reverse to match the reversed factor levels
names(watershed_colors) <- rev(levels(boxplot_data$mgmt_river))

# Create sleek, modern boxplot
p <- ggplot(boxplot_data, aes(x = mgmt_river, y = iso_pred, fill = mgmt_river)) +
  geom_boxplot(
    alpha = 0.85,
    outlier.size = 2,
    outlier.shape = 21,
    outlier.alpha = 0.6,
    outlier.fill = "gray40",
    color = "gray30",
    linewidth = 0.4,
    width = 0.7
  ) +
  scale_fill_manual(values = watershed_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0.704, 0.713), 
    breaks = seq(0.704, 0.713, 0.001),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_flip() +
  labs(
    title = "Isotope Ratio Distribution by Management Unit",
    subtitle = "Red = Upstream • Blue = Downstream | Only streams used in assignment (with all priors applied)",
    x = NULL,
    y = "iso_pred (δ²H)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0, margin = margin(b = 10)),
    axis.text.y = element_text(size = 11, color = "gray20"),
    axis.text.x = element_text(size = 10, color = "gray30"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# Save plot
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(OUTPUT_DIR, "mgmt_unit_iso_pred_boxplot.png")
ggsave(output_file, p, width = 10, height = 8, dpi = 300)
cat(paste("Saved plot to:", output_file, "\n"))

################################################################################
# EXPORT SUMMARY TABLE
################################################################################

# Save summary table
csv_file <- file.path(OUTPUT_DIR, "mgmt_unit_iso_pred_summary.csv")
write.csv(mgmt_iso_summary, csv_file, row.names = FALSE)
cat(paste("Saved summary table to:", csv_file, "\n"))

# Print summary
cat("\n=== SUMMARY ===\n")
print(mgmt_iso_summary, n = Inf)
cat("\n")