# HUC Protection Comparison Plot (30% Baseline)
# This script creates a single plot showing protection amount over time
# for all HUCs with a 30% baseline harvest rate

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Read the HUC closure protection data
huc_protection <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/huc_closure_summary_Kusko_2017_to_2021_06-01_to_06-11.csv")

# Rename column if needed
huc_protection <- huc_protection %>%
  rename(HUC.Protection = HUC.Protection..)

# Set the baseline harvest rate to 30%
baseline_rate <- 0.3  # 30%

# Calculate protection amount for each HUC and year
protection_data <- huc_protection %>%
  mutate(
    Baseline_Harvest_Rate = baseline_rate,
    Total_Harvest_Rate = (1 - HUC.Protection/100) * baseline_rate,
    Protection_Amount = baseline_rate - Total_Harvest_Rate
  )

# Create output directory for visualizations
dir.create("HUC_Protection_Visualizations", showWarnings = FALSE)

# Save the calculated data for reference
write.csv(protection_data, "HUC_Protection_Visualizations/Protection_Data_30pct_Baseline.csv", row.names = FALSE)

# Generate a color palette with enough distinct colors for all HUCs
n_hucs <- length(unique(protection_data$HUC.Name))
huc_colors <- brewer.pal(min(n_hucs, 8), "Dark2")  # Use Dark2 palette for better distinction
# If we have more than 8 HUCs, extend the palette
if (n_hucs > 8) {
  huc_colors <- colorRampPalette(huc_colors)(n_hucs)
}

# Create the plot
p <- ggplot(protection_data, aes(x = Year, y = Protection_Amount, 
                                 group = HUC.Name, color = HUC.Name)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3.5) +
  scale_color_manual(values = huc_colors, name = "HUC Name") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Harvest Rate Reduction Due to Protection",
       subtitle = paste0("Comparison across all HUCs with ", baseline_rate*100, "% baseline harvest rate"),
       x = "Year",
       y = "Amount of Protection (Harvest Rate Reduction)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave("HUC_Protection_Visualizations/Protection_Comparison_30pct_Baseline.png", 
       p, width = 12, height = 8, dpi = 300)

