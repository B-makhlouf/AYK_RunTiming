# Portfolio Effect Figure — animated GIF
# Requires: ggplot2, gganimate, gifski, dplyr
# Install once with:
#   install.packages(c("ggplot2", "gganimate", "gifski", "dplyr"))

library(ggplot2)
library(gganimate)
library(dplyr)

set.seed(42)

pop_colors <- c(
  "#2166AC", "#4DD8D0", "#CC0000", "#E87722", "#6A0DAD", "#2CA02C",
  "#AACC00", "#E040FB", "#8B4513", "#1E90FF", "#A67C00", "#FF6699"
)

sds  <- c(0.18, 0.06, 0.22, 0.12, 0.08, 0.25, 0.10, 0.20, 0.15, 0.05, 0.17, 0.09)
mean_val <- 0.5
years    <- 1:10

# Build data frame
df <- do.call(rbind, lapply(seq_along(pop_colors), function(i) {
  vals <- mean_val + rnorm(length(years), sd = sds[i])
  vals <- pmax(pmin(vals, 0.98), 0.02)
  data.frame(
    pop   = factor(i),
    year  = years,
    abund = vals,
    color = pop_colors[i]
  )
}))

# The animation reveals one year at a time; transition_reveal draws lines
p <- ggplot(df, aes(x = year, y = abund, group = pop, color = pop)) +
  geom_hline(yintercept = mean_val, linetype = "dashed",
             colour = "grey50", linewidth = 0.6, alpha = 0.5) +
  geom_line(linewidth = 1.6, alpha = 0.85) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = setNames(pop_colors, levels(df$pop))) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    panel.grid      = element_blank()
  ) +
  transition_reveal(year)   # <-- key gganimate call

# Render and save
anim <- animate(
  p,
  nframes  = 60,
  fps      = 6,
  width    = 900,
  height   = 550,
  renderer = gifski_renderer("portfolio_effect_animated.gif")
)

anim_save("portfolio_effect_animated.gif", anim)
message("Done! GIF saved to working directory.")
