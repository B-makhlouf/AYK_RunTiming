# Simple Run Duration Analysis: 5% to 100% by Year
# Shows duration values within each single year, organized by watershed position

# Load required libraries
library(tidyverse)
library(ggplot2)
library(viridis)

# Read the data
data_path <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Cumulative_Distribution/CSV/cumulative_distribution_wide_Kusko_3day_intervals.csv"
df <- read.csv(data_path)

# Define watershed ordering (upstream to downstream)
get_watershed_order <- function() {
  watershed_order <- c(
    "N..Fork.Kusko",
    "E..Fork.Kuskokwim.River", 
    "S..Fork.Kusko",
    "Upper.Kusko.Main",
    "Big.River",
    "Takotna.and.Nixon.Fork",
    "Tatlawiksuk",
    "Swift",
    "Stony", 
    "Holitna",
    "Hoholitna",
    "Middle.Kusko.Main",
    "George",
    "Oskakawlik", 
    "Holokuk",
    "Aniak",
    "Tuluksak",
    "Kisaralik",
    "Kwethluk",
    "Johnson",
    "Lower.Kusko"
  )
  return(watershed_order)
}

# Function to interpolate day of year for a given percentage
interpolate_doy <- function(cumulative_values, doy_values, target_percentage) {
  valid_idx <- !is.na(cumulative_values)
  if(sum(valid_idx) < 2) return(NA)
  
  cumulative_clean <- cumulative_values[valid_idx]
  doy_clean <- doy_values[valid_idx]
  
  if(cumulative_clean[1] >= target_percentage) {
    return(doy_clean[1])
  }
  
  if(max(cumulative_clean) < target_percentage) {
    return(NA)
  }
  
  idx <- which(cumulative_clean >= target_percentage)[1]
  if(idx == 1) {
    return(doy_clean[1])
  }
  
  # Linear interpolation
  x1 <- cumulative_clean[idx-1]
  x2 <- cumulative_clean[idx]
  y1 <- doy_clean[idx-1]
  y2 <- doy_clean[idx]
  
  interpolated_doy <- y1 + (target_percentage - x1) * (y2 - y1) / (x2 - x1)
  return(interpolated_doy)
}

# Calculate 5% and 100% days for each stock and year
stock_columns <- grep("^mgmt_", names(df), value = TRUE)
results <- data.frame()

for(stock in stock_columns) {
  stock_name <- gsub("^mgmt_", "", stock)
  
  for(year in unique(df$year)) {
    year_data <- df[df$year == year, ]
    cumulative_values <- year_data[[stock]]
    doy_values <- year_data$doy
    
    # Calculate 5% and 100% days
    day_5 <- interpolate_doy(cumulative_values, doy_values, 1)
    day_100 <- interpolate_doy(cumulative_values, doy_values, 50)
    
    # Calculate duration (days between 5% and 100%)
    duration <- if(!is.na(day_5) && !is.na(day_100)) day_100 - day_5 else NA
    
    results <- rbind(results, data.frame(
      year = year,
      stock = stock_name,
      duration_days = duration
    ))
  }
}

# Remove rows with missing duration data
results_clean <- results[!is.na(results$duration_days), ]

# Clean stock names and apply watershed ordering
results_clean$stock_clean <- gsub("\\.", " ", results_clean$stock)

# Handle naming variations
name_mapping <- c(
  "N  Fork Kusko" = "N. Fork Kusko",
  "E  Fork Kuskokwim River" = "E. Fork Kuskokwim River",
  "S  Fork Kusko" = "S. Fork Kusko"
)

for(old_name in names(name_mapping)) {
  results_clean$stock_clean[results_clean$stock_clean == old_name] <- name_mapping[old_name]
}

# Apply watershed ordering
watershed_order <- get_watershed_order()
watershed_order_clean <- gsub("\\.", " ", watershed_order)

for(old_name in names(name_mapping)) {
  watershed_order_clean[watershed_order_clean == old_name] <- name_mapping[old_name]
}

stocks_in_data <- unique(results_clean$stock_clean)
final_order <- watershed_order_clean[watershed_order_clean %in% stocks_in_data]

# Add any missing stocks
missing_stocks <- setdiff(stocks_in_data, final_order)
if(length(missing_stocks) > 0) {
  final_order <- c(final_order, missing_stocks)
}

# Apply ordering (reverse for coord_flip so upstream is at top)
results_clean$stock_ordered <- factor(results_clean$stock_clean, levels = rev(final_order))

# CREATE SIMPLE PLOT: One panel per year, watershed ordered
p_simple <- ggplot(results_clean, aes(x = stock_ordered, y = duration_days)) +
  geom_col(aes(fill = duration_days), alpha = 0.8, color = "white", size = 0.3) +
  facet_wrap(~year, ncol = 3) +
  scale_fill_viridis_c(name = "Days", option = "plasma") +
  coord_flip() +
  labs(
    title = "Run Duration (5% to 100%) by Management Unit",
    subtitle = "Days from start to end of run | Ordered by watershed position (upstream → downstream)",
    x = "Management Unit",
    y = "Duration (Days)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

# Print and save
print(p_simple)
ggsave("run_duration_simple_by_year.png", p_simple, 
       width = 14, height = 10, dpi = 300, bg = "white")

# Print the data for verification
cat("Sample of calculated durations:\n")
print(head(results_clean, 10))

cat("\nWatershed order applied (upstream to downstream):\n")
for(i in 1:length(final_order)) {
  cat(paste(i, ".", final_order[i], "\n"))
}

cat("\nSimple plot saved as: run_duration_simple_by_year.png\n")