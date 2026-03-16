# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)

# Read the data
data <- read_excel("/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/AYKEscapement.xlsx")

# View the structure of the data
str(data)
head(data)

# ===== ESCAPEMENT FIGURE =====

# Create a 3-panel figure showing escapement per year for each watershed
p_escapement <- ggplot(data, aes(x = Year, y = Escapement)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_point(size = 2, color = "darkgreen", alpha = 0.7) +
  facet_wrap(~ River, ncol = 1, scales = "free_y") +
  labs(
    title = "Chinook Salmon Escapement by Watershed",
    x = "Year",
    y = "Escapement (number of fish)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = scales::comma)

# Display the plot
print(p_escapement)

# Save the plot
ggsave(
  filename = "/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/salmon_escapement_3panel.png",
  plot = p_escapement,
  width = 10,
  height = 8,
  dpi = 300
)

# ===== HARVEST FIGURE =====

# Create a 3-panel figure showing harvest per year for each watershed
p_harvest <- ggplot(data, aes(x = Year, y = Harvest)) +
  geom_line(linewidth = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick", alpha = 0.7) +
  facet_wrap(~ River, ncol = 1, scales = "free_y") +
  labs(
    title = "Chinook Salmon Harvest by Watershed",
    x = "Year",
    y = "Harvest (number of fish)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = scales::comma)

# Display the plot
print(p_harvest)

# Save the plot
ggsave(
  filename = "/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/salmon_harvest_3panel.png",
  plot = p_harvest,
  width = 10,
  height = 8,
  dpi = 300
)

# ===== COMPARISON PLOTS (BONUS) =====

# Escapement comparison - all three watersheds on same panel
p_escapement_comparison <- ggplot(data, aes(x = Year, y = Escapement, color = River)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "Chinook Salmon Escapement Comparison",
    x = "Year",
    y = "Escapement (number of fish)",
    color = "Watershed"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_brewer(palette = "Set1")

print(p_escapement_comparison)

ggsave(
  filename = "/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/salmon_escapement_comparison.png",
  plot = p_escapement_comparison,
  width = 10,
  height = 6,
  dpi = 300
)

# Harvest comparison - all three watersheds on same panel
p_harvest_comparison <- ggplot(data, aes(x = Year, y = Harvest, color = River)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "Chinook Salmon Harvest Comparison",
    x = "Year",
    y = "Harvest (number of fish)",
    color = "Watershed"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_brewer(palette = "Set1")

print(p_harvest_comparison)

ggsave(
  filename = "/Users/benjaminmakhlouf/Research_repos/Schindler_GitHub/Arctic_Yukon_Kuskokwim_Data/salmon_harvest_comparison.png",
  plot = p_harvest_comparison,
  width = 10,
  height = 6,
  dpi = 300
)

# ===== SUMMARY STATISTICS =====

# Print summary statistics by watershed for Escapement
cat("\n=== Escapement Summary Statistics by Watershed ===\n\n")
escapement_stats <- data %>%
  group_by(River) %>%
  summarize(
    Years = n(),
    Min_Escapement = min(Escapement, na.rm = TRUE),
    Max_Escapement = max(Escapement, na.rm = TRUE),
    Mean_Escapement = mean(Escapement, na.rm = TRUE),
    Median_Escapement = median(Escapement, na.rm = TRUE),
    SD_Escapement = sd(Escapement, na.rm = TRUE)
  )

print(escapement_stats)

# Print summary statistics by watershed for Harvest
cat("\n=== Harvest Summary Statistics by Watershed ===\n\n")
harvest_stats <- data %>%
  group_by(River) %>%
  summarize(
    Years = n(),
    Min_Harvest = min(Harvest, na.rm = TRUE),
    Max_Harvest = max(Harvest, na.rm = TRUE),
    Mean_Harvest = mean(Harvest, na.rm = TRUE),
    Median_Harvest = median(Harvest, na.rm = TRUE),
    SD_Harvest = sd(Harvest, na.rm = TRUE)
  )

print(harvest_stats)

# Calculate and display harvest rate (harvest/total return) by watershed
cat("\n=== Average Harvest Rate by Watershed ===\n\n")
harvest_rate_stats <- data %>%
  mutate(Harvest_Rate = Harvest / Total_Return * 100) %>%
  group_by(River) %>%
  summarize(
    Mean_Harvest_Rate_Percent = mean(Harvest_Rate, na.rm = TRUE),
    SD_Harvest_Rate = sd(Harvest_Rate, na.rm = TRUE),
    Min_Harvest_Rate = min(Harvest_Rate, na.rm = TRUE),
    Max_Harvest_Rate = max(Harvest_Rate, na.rm = TRUE)
  )

print(harvest_rate_stats)

cat("\n=== Plots saved successfully! ===\n")