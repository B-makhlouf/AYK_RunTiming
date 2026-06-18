library(ggplot2)
library(dplyr)

set.seed(42)

# ── Shared parameters ──────────────────────────────────────────────────────────

pop_colors <- c(
  "#2166AC", "#4DD8D0", "#CC0000", "#E87722", "#6A0DAD", "#2CA02C",
  "#AACC00", "#E040FB", "#8B4513", "#1E90FF", "#A67C00", "#FF6699"
)

sds      <- c(0.18, 0.06, 0.22, 0.12, 0.08, 0.25, 0.10, 0.20, 0.15, 0.05, 0.17, 0.09)
mean_val <- 0.5
years    <- 1:10
shared   <- rnorm(length(years), sd = 0.12)

# ── Build data ─────────────────────────────────────────────────────────────────

# Individual population series
pop_df <- do.call(rbind, lapply(seq_along(pop_colors), function(i) {
  vals <- mean_val + shared + rnorm(length(years), sd = sds[i])
  vals <- pmax(pmin(vals, 0.98), 0.02)
  data.frame(pop = factor(i), year = years, abundance = vals)
}))

# Aggregate (mean across populations each year)
agg_df <- pop_df |>
  group_by(year) |>
  summarise(abundance = mean(abundance), .groups = "drop")

# ── Shared theme ───────────────────────────────────────────────────────────────

base_theme <- theme_grey(base_size = 13) +
  theme(
    legend.position  = "none",
    axis.text        = element_blank(),
    axis.ticks       = element_blank()
  )

y_scale <- scale_y_continuous(limits = c(0, 1.05),
                               breaks = c(0, 0.25, 0.5, 0.75, 1.0))

# ── Figure 1: individual populations ──────────────────────────────────────────

p1 <- ggplot(pop_df, aes(x = year, y = abundance, group = pop, colour = pop)) +
  geom_hline(yintercept = mean_val, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_line(linewidth = 1.4, alpha = 0.85) +
  scale_colour_manual(values = setNames(pop_colors, levels(pop_df$pop))) +
  scale_x_continuous(breaks = years) +
  y_scale +
  labs(x = "Year", y = "Relative abundance") +
  base_theme

ggsave("C:/Users/makhl/Research Repos/AYK_RunTiming/Figures/PresFigures/portfolio_populations.png",
       plot = p1, width = 8, height = 5, dpi = 200, bg = "white")

# ── Figure 2: aggregate only ───────────────────────────────────────────────────

p2 <- ggplot(agg_df, aes(x = year, y = abundance)) +
  geom_hline(yintercept = mean_val, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_line(colour = "#C97B7E", linewidth = 2.0, alpha = 0.9) +
  scale_x_continuous(breaks = years) +
  y_scale +
  labs(x = "Year", y = "Relative abundance") +
  base_theme

ggsave("C:/Users/makhl/Research Repos/AYK_RunTiming/Figures/PresFigures/portfolio_aggregate.png",
       plot = p2, width = 8, height = 5, dpi = 200, bg = "white")

message("Saved to C:/Users/makhl/Research Repos/AYK_RunTiming/Figures/PresFigures/")
