################################################################################
# FRONT-END CLOSURE BOXPLOT - CORRECT CALCULATION WITH WATERSHED ORDERING
################################################################################
# PURPOSE: Shows what % of EACH management unit's total annual production 
#          falls within the Q1 closure window (start through June 11)
# NOTE: Now uses standardized watershed ordering (upstream to downstream)
################################################################################

# Load libraries
library(dplyr)
library(ggplot2)

# Load required functions - use relative path or check if functions exist
if(!exists("get_watershed_order")) {
  # Define watershed ordering functions directly in this script
  get_watershed_order <- function() {
    c("N. Fork Kusko", "E. Fork Kuskokwim River", "S. Fork Kusko",
      "Upper Kusko Main", "Big River", "Takotna and Nixon Fork",
      "Tatlawiksuk", "Swift", "Stony", "Holitna and Hoholitna",
      "Middle Kusko Main", "George", "Oskakawlik", "Holokuk",
      "Aniak", "Tuluksak", "Kisaralik", "Kwethluk", "Johnson", "Lower Kusko")
  }
  
  apply_watershed_order <- function(data, mgmt_col = "mgmt_river", reverse_for_plots = FALSE) {
    standard_order <- get_watershed_order()
    if (reverse_for_plots) standard_order <- rev(standard_order)
    
    units_in_data <- unique(data[[mgmt_col]])
    final_order <- standard_order[standard_order %in% units_in_data]
    
    missing_units <- setdiff(units_in_data, standard_order)
    if (length(missing_units) > 0) {
      warning("Found management units not in standard order: ", paste(missing_units, collapse = ", "))
      final_order <- c(final_order, missing_units)
    }
    
    data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
    return(data)
  }
}

# Output directory
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front_End_Closure_Boxplots"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# STEP 1: LOAD ALL QUARTILE DATA
################################################################################

# Load ALL quartile data (Q1, Q2, Q3, Q4) to calculate totals
mgmt_data <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")

# Remove Johnson river and clean quartile names
mgmt_data_clean <- mgmt_data %>%
  filter(mgmt_river != "Johnson") %>%
  mutate(quartile_clean = case_when(
    quartile == "Q1" ~ "Q1",
    quartile == "Q2" ~ "Q2", 
    quartile == "Q3" ~ "Q3",
    quartile == "Q4" ~ "Q4",
    TRUE ~ quartile
  ))

print(paste("Loaded data for", length(unique(mgmt_data_clean$mgmt_river)), "management units"))
print(paste("Years available:", paste(sort(unique(mgmt_data_clean$year)), collapse = ", ")))
print(paste("Quartiles available:", paste(sort(unique(mgmt_data_clean$quartile_clean)), collapse = ", ")))

################################################################################
# STEP 2: CALCULATE TOTAL ANNUAL PRODUCTION PER UNIT
################################################################################

# Calculate total annual production for each management unit in each year
# Sum the total_run_prop across all quartiles (Q1+Q2+Q3+Q4)
annual_totals <- mgmt_data_clean %>%
  group_by(year, mgmt_river) %>%
  summarise(
    total_annual_production = sum(total_run_prop, na.rm = TRUE),
    .groups = "drop"
  )

print("\nSample annual totals:")
print(head(annual_totals))

################################################################################
# STEP 3: CALCULATE Q1 CLOSURE PERCENTAGE
################################################################################

# Get Q1 production values
q1_production <- mgmt_data_clean %>%
  filter(quartile_clean == "Q1") %>%
  select(year, mgmt_river, q1_production = total_run_prop)

# Join Q1 with annual totals and calculate closure percentage
boxplot_data <- q1_production %>%
  left_join(annual_totals, by = c("year", "mgmt_river")) %>%
  mutate(
    closure_percentage = (q1_production / total_annual_production) * 100
  ) %>%
  select(year, mgmt_river, q1_production, total_annual_production, closure_percentage)

# Verify calculation with examples
print("\n=== CALCULATION VERIFICATION ===")
print("Sample calculations showing the formula:")
sample_calcs <- boxplot_data %>% 
  slice_head(n = 5) %>%
  mutate(
    formula_check = round((q1_production / total_annual_production) * 100, 2)
  )
print(sample_calcs)

# Check for any issues
missing_data <- boxplot_data %>% filter(is.na(closure_percentage))
if(nrow(missing_data) > 0) {
  print("\nWarning: Missing data found:")
  print(missing_data)
}

################################################################################
# STEP 4: APPLY WATERSHED ORDERING
################################################################################

# Apply watershed ordering (reverse for coord_flip so upstream appears at top)
boxplot_data <- apply_watershed_order(boxplot_data, "mgmt_river", reverse_for_plots = TRUE)

# Print the final ordering
print("\nWatershed order applied (will appear top to bottom in plot):")
final_levels <- levels(boxplot_data$mgmt_river)
for (i in 1:length(final_levels)) {
  cat(paste(i, ".", final_levels[i], "\n"))
}

################################################################################
# STEP 5: CREATE BOXPLOT
################################################################################

year_range <- paste0(min(boxplot_data$year), "-", max(boxplot_data$year))

boxplot <- ggplot(boxplot_data, aes(x = mgmt_river, y = closure_percentage)) +
  geom_boxplot(
    fill = "lightblue", 
    alpha = 0.7, 
    outlier.size = 2.5,
    outlier.shape = 16,
    linewidth = 0.6,
    color = "darkblue"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, max(boxplot_data$closure_percentage, na.rm = TRUE) * 1.05),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "Front-End Closure Protection by Management Unit",
    subtitle = paste("% of EACH UNIT'S total annual production within Q1 closure window | Watershed order: upstream → downstream | Years:", year_range),
    x = "Management Unit (Watershed Position: Upstream → Downstream)",
    y = "% of Unit's Total Annual Production in Q1 Closure Window",
    caption = "Shows what % of each unit's total annual production (Q1+Q2+Q3+Q4) occurs within Q1 closure period\nOrdered by position in watershed from headwaters (top) to mouth (bottom)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey20"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray50"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

################################################################################
# STEP 6: SAVE PLOT AND SUMMARY
################################################################################

ggsave(file.path(OUTPUT_DIR, "front_end_closure_protection_boxplot_WATERSHED_ORDERED.png"), 
       boxplot, width = 12, height = 10, dpi = 300, bg = "white")

print(boxplot)

# Summary statistics with watershed ordering maintained
summary_stats <- boxplot_data %>%
  group_by(mgmt_river) %>%
  summarise(
    mean_closure_pct = round(mean(closure_percentage, na.rm = TRUE), 1),
    median_closure_pct = round(median(closure_percentage, na.rm = TRUE), 1),
    min_closure_pct = round(min(closure_percentage, na.rm = TRUE), 1),
    max_closure_pct = round(max(closure_percentage, na.rm = TRUE), 1),
    n_years = n(),
    .groups = "drop"
  ) %>%
  # Maintain the watershed ordering in the summary
  arrange(mgmt_river)

write.csv(summary_stats, file.path(OUTPUT_DIR, "closure_protection_summary_WATERSHED_ORDERED.csv"), row.names = FALSE)

print("\n=== SUMMARY STATISTICS (WATERSHED ORDERED) ===")
print("Top 5 management units (upstream first):")
print(head(summary_stats))

print("\n=== ANALYSIS COMPLETE ===")
print("✓ Applied standardized watershed ordering (upstream → downstream)")
print("✓ Upstream management units appear at TOP of plot")
print("✓ Downstream management units appear at BOTTOM of plot")
print(paste("✓ Plot saved to:", file.path(OUTPUT_DIR, "front_end_closure_protection_boxplot_WATERSHED_ORDERED.png")))