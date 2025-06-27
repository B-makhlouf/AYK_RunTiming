################################################################################
# FRONT-END CLOSURE BOXPLOT - CORRECT CALCULATION
################################################################################
# PURPOSE: Shows what % of EACH management unit's total annual production 
#          falls within the Q1 closure window (start through June 11)
#
# CORRECT CALCULATION:
# For each management unit in each year:
# 1. Get Q1 production value (total_run_prop for Q1)
# 2. Get total annual production (sum of Q1+Q2+Q3+Q4 total_run_prop)
# 3. closure_percentage = (Q1 production / Total annual production) × 100
#
# INTERPRETATION:
# - High % = Most of that unit's fish arrive during Q1 closure (early timing)
# - Low % = Most of that unit's fish arrive after Q1 closure (late timing)
################################################################################

# Load libraries
library(dplyr)
library(ggplot2)

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
# STEP 4: SET WATERSHED ORDER (UPSTREAM TO DOWNSTREAM)
################################################################################

# Define watershed position order (same as average production analysis)
watershed_order <- c(
  "N. Fork Kusko", "E. Fork Kuskokwim", "S. Fork Kusko", 
  "Takotna and Nixon Fork", "Big River", "Upper Kusko Main",
  "Tatlawiksuk", "Kwethluk", "Stony", "Swift",
  "Holitna and Hoholitna", "George", "Oskakawlik", 
  "Middle Kusko Main", "Holokuk", "Aniak", "Tuluksak",
  "Kisaralik", "Hoholitna", "Johnson", "Lower Kusko"
)

# Filter to units present in data, maintain watershed order
units_in_data <- unique(boxplot_data$mgmt_river)
final_order <- watershed_order[watershed_order %in% units_in_data]

# Apply ordering (reverse for coord_flip so upstream appears at top)
boxplot_data$mgmt_river <- factor(boxplot_data$mgmt_river, levels = rev(final_order))

print(paste("\nWatershed order (upstream to downstream):", paste(final_order, collapse = " → ")))

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
    x = "Management Unit (Watershed Position)",
    y = "% of Unit's Total Annual Production in Q1 Closure Window",
    caption = "Shows what % of each unit's total annual production (Q1+Q2+Q3+Q4) occurs within Q1 closure period"
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

ggsave(file.path(OUTPUT_DIR, "front_end_closure_protection_boxplot_CORRECTED.png"), 
       boxplot, width = 12, height = 10, dpi = 300, bg = "white")

print(boxplot)

# Summary statistics
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
  mutate(mgmt_river = factor(mgmt_river, levels = final_order)) %>%
  arrange(mgmt_river)

write.csv(summary_stats, file.path(OUTPUT_DIR, "closure_protection_summary_CORRECTED.csv"), row.names = FALSE)


