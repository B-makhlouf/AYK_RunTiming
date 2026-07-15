################################################################################
# SUPPLEMENTAL FIGURE: DAILY CPUE AND SELECTED SAMPLES BY YEAR
################################################################################
# PURPOSE: For each year, show day-of-year (DOY), daily CPUE proportion, and the
#          number of selected otolith samples per day, in a single faceted figure.
# OUTPUT:  - supplemental_cpue_samples_by_year.csv (doy, cpue, n_samples, year)
#          - FigureS_CPUE_samples_by_year.png (6 panels, one per year)
################################################################################

library(dplyr)
library(readr)
library(ggplot2)

if (!exists("PROJECT_ROOT")) source("config.R")

years     <- CFG$years
watershed <- "Kusko"
output_dir <- CPUE_DIR
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# BUILD PER-DAY TABLE
################################################################################
# For each day: daily CPUE proportion is constant within a DOY (take first);
# selected samples = number of fish records (non-missing Fish_id) on that day.

all_years <- list()
for (year in years) {
  file_path <- file.path(
    NATAL_DATA_DIR,
    paste0(year, "_", watershed, "_Natal_Origins_Genetics_CPUE.csv")
  )
  if (!file.exists(file_path)) { warning("Missing file for ", year); next }

  daily <- read.csv(file_path, stringsAsFactors = FALSE) %>%
    group_by(doy = DOY) %>%
    summarise(
      cpue      = dplyr::first(dailyCPUEprop),
      n_samples = sum(!is.na(Fish_id)),
      .groups   = "drop"
    ) %>%
    mutate(cpue = ifelse(is.na(cpue), 0, cpue), year = year) %>%
    arrange(doy)

  all_years[[as.character(year)]] <- daily
}

data <- bind_rows(all_years) %>%
  select(year, doy, cpue, n_samples)

write_csv(data, file.path(output_dir, "supplemental_cpue_samples_by_year.csv"))

################################################################################
# FIGURE: bars = selected samples (right axis), line = daily CPUE (left axis)
################################################################################

line_col <- "#2c6fbb"
bar_col  <- "#c9ced6"

# Scale sample counts onto the CPUE axis for dual-axis display
scale_factor <- max(data$cpue) / max(data$n_samples)

p <- ggplot(data, aes(x = doy)) +
  geom_col(aes(y = n_samples * scale_factor),
           fill = bar_col, width = 0.9) +
  geom_line(aes(y = cpue), color = line_col, linewidth = 0.7) +
  facet_wrap(~ year, ncol = 1, scales = "fixed") +
  scale_y_continuous(
    name = "Daily CPUE proportion",
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of selected samples")
  ) +
  labs(x = "Day of year (DOY)",
       title = "Daily CPUE and selected otolith samples by day of year, 2017-2022") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text        = element_text(face = "bold", hjust = 0),
    axis.title.y      = element_text(color = line_col),
    axis.text.y.left  = element_text(color = line_col),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(output_dir, "FigureS_CPUE_samples_by_year.png"),
       p, width = 8.5, height = 11, dpi = 300, bg = "white")

cat("Saved supplemental CSV and figure to", output_dir, "\n")
