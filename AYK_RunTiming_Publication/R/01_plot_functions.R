################################################################################
# 05_VISUALIZATION_FUNCTIONS.R - ALL PLOTTING AND MAPPING FUNCTIONS
################################################################################
# Shared visualization functions used across all analyses
# Load this after 00_setup.R
################################################################################

cat("Loading visualization functions...\n")

################################################################################
# TRIBUTARY MAPPING FUNCTIONS
################################################################################

#' Create DOY histogram for tributary maps
create_doy_histogram <- function(full_dataset, current_subset, title = NULL) {
  if (is.null(title)) title <- "DOY Distribution"
  
  highlight_color <- "tomato"
  background_color <- "gray70"
  
  doy_to_date <- function(doy, year = 2024) {
    as.Date(doy - 1, origin = paste0(year, "-01-01"))
  }
  
  doy_breaks <- seq(140, 210, by = 10)
  
  ggplot() + 
    # Background curve
    geom_line(data = full_dataset, aes(x = DOY, y = dailyCPUEprop), 
              color = "gray40", linewidth = 1, alpha = 0.5) +
    geom_ribbon(data = full_dataset, aes(x = DOY, ymin = 0, ymax = dailyCPUEprop), 
                fill = background_color, alpha = 0.3) +
    
    # Highlighted subset
    geom_line(data = current_subset, aes(x = DOY, y = dailyCPUEprop), 
              color = "black", linewidth = 2) +
    geom_ribbon(data = current_subset, aes(x = DOY, ymin = 0, ymax = dailyCPUEprop), 
                fill = highlight_color, alpha = 0.7) +
    
    # Boundary lines
    geom_vline(xintercept = c(min(current_subset$DOY), max(current_subset$DOY)), 
               linetype = "dashed", color = "darkred") +
    
    scale_x_continuous(
      limits = c(140, 210),
      breaks = doy_breaks,
      labels = function(x) paste0(x, "\n", format(doy_to_date(x), "%b %d"))
    ) +
    scale_y_continuous(limits = c(0, 0.1)) +
    coord_cartesian(xlim = c(140, 210), ylim = c(0, 0.1), expand = FALSE) +
    
    labs(
      title = title,
      subtitle = "Current subset highlighted in red",
      x = "Day of Year (Date)", 
      y = "Daily CPUE Proportion"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

#' Create tributary map (matching original style exactly)
create_tributary_map <- function(basin, edges, basin_assign_norm, year, watershed, 
                                 subset_label, output_path, StreamOrderPrior = NULL, pid_prior = NULL, gg_hist = NULL) {
  
  png(file = output_path, width = 9, height = 8, units = "in", res = 300, bg = "white")
  
  # Use YlOrRd palette with 9 colors expanded to 10 (matching original)
  pallete <- brewer.pal(9, "YlOrRd")
  pallete_expanded <- colorRampPalette(pallete)(10)
  
  # Color coding with bins at every 0.1 (matching original exactly)
  colcode <- rep("gray60", length(basin_assign_norm))
  colcode[basin_assign_norm == 0] <- 'white'
  colcode[basin_assign_norm > 0 & basin_assign_norm <= 0.2] <- pallete_expanded[1]
  colcode[basin_assign_norm > 0.2 & basin_assign_norm <= 0.4] <- pallete_expanded[4]
  colcode[basin_assign_norm > 0.4 & basin_assign_norm <= 0.6] <- pallete_expanded[5]
  colcode[basin_assign_norm > 0.6 & basin_assign_norm <= 0.7] <- pallete_expanded[7]
  colcode[basin_assign_norm > 0.7 & basin_assign_norm <= 0.8] <- pallete_expanded[8]
  colcode[basin_assign_norm > 0.8 & basin_assign_norm <= 0.9] <- pallete_expanded[9]
  colcode[basin_assign_norm > 0.9 & basin_assign_norm <= 1.0] <- pallete_expanded[10]
  
  # Override with grays for excluded areas (if priors provided)
  if (!is.null(StreamOrderPrior)) {
    colcode[which(StreamOrderPrior == 0)] <- 'gray60'
  }
  if (!is.null(pid_prior)) {
    colcode[which(pid_prior == 0)] <- 'gray60'
  }
  
  # Set linewidths based on stream order (matching original exactly)
  stream_order_lwd <- edges$Str_Order
  linewidths <- rep(1, length(stream_order_lwd))
  
  if (watershed == "Yukon") {
    linewidths <- ifelse(stream_order_lwd == 9, 3.7, linewidths)
    linewidths <- ifelse(stream_order_lwd == 8, 2.5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 7, 1.7, linewidths)
    linewidths <- ifelse(stream_order_lwd == 6, 1.5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 5, 1, linewidths)
    linewidths <- ifelse(stream_order_lwd == 4, 1, linewidths)
    linewidths <- ifelse(stream_order_lwd == 3, 1, linewidths)
  } else {
    linewidths <- ifelse(stream_order_lwd == 9, 5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 8, 4, linewidths)
    linewidths <- ifelse(stream_order_lwd == 7, 3, linewidths)
    linewidths <- ifelse(stream_order_lwd == 6, 2, linewidths)
    linewidths <- ifelse(stream_order_lwd == 5, 1.8, linewidths)
    linewidths <- ifelse(stream_order_lwd == 4, 1.5, linewidths)
    linewidths <- ifelse(stream_order_lwd == 3, 1, linewidths)
  }
  
  # Generate title (matching original format)
  plot_title <- paste0(subset_label, "\nYear:", year, " River:", watershed)
  
  # Set plot margins (matching original)
  par(mar = c(8, 4, 4, 2), bg = "white")
  
  # Plot basin and edges (matching original exactly)
  plot(st_geometry(basin), col = 'gray60', border = 'gray60', main = plot_title, bg = "white")
  plot(st_geometry(edges), col = colcode, pch = 16, axes = FALSE, add = TRUE, lwd = linewidths)
  
  # Add legend (matching original exactly)
  legend("topleft", 
         legend = c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", 
                    "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"), 
         col = pallete_expanded, 
         lwd = 5, 
         title = "Relative posterior density", 
         bty = "n",
         bg = "white")
  
  # ADD HISTOGRAM OVERLAY (THIS WAS MISSING) - EXACTLY LIKE ORIGINAL
  if (!is.null(gg_hist)) {
    # Modify the histogram specifically for grid viewport use (from original code)
    limited_hist <- gg_hist +
      scale_x_continuous(limits = c(140, 200)) +
      scale_y_continuous(limits = c(0, 0.1)) +
      coord_cartesian(xlim = c(140, 200), ylim = c(0, 0.1), expand = FALSE) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(0, 0, 0, 0)
      )
    
    # Create viewport with explicit scaling (from original code)
    vp_hist <- viewport(
      x = 0.5, y = 0.05, 
      width = 0.7, height = 0.2, 
      just = c("center", "bottom")
    )
    
    # Print the modified histogram (from original code)
    print(limited_hist, vp = vp_hist)
  }
  
  dev.off()
  
  # Reset par to default
  par(mar = c(5, 4, 4, 2) + 0.1, bg = "white")
  
  message(paste("Created tributary map:", basename(output_path)))
}

#' Create management river map (for DOY quartile analysis)
create_mgmt_river_map <- function(mgmt_result, edges, basin, year, watershed, 
                                  subset_label, output_path, gg_hist = NULL) {
  
  png(file = output_path, width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Set up plotting layout  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2, 
                                             heights = unit(c(0.7, 0.3), "npc"),
                                             widths = unit(c(0.6, 0.4), "npc"))))
  
  # Add production proportion to edges
  edges$production_proportion <- NA
  for (i in 1:nrow(mgmt_result)) {
    mgmt_name <- mgmt_result$mgmt_river[i]
    prod_prop <- mgmt_result$production_proportion[i]
    edges$production_proportion[edges$mgmt_river == mgmt_name] <- prod_prop
  }
  
  # Filter to only managed edges
  managed_edges <- edges[!is.na(edges$mgmt_river) & edges$mgmt_river != "", ]
  
  # Convert to sf objects for ggplot
  managed_edges_sf <- st_as_sf(managed_edges)
  basin_sf <- st_as_sf(basin)
  
  # Main map plot
  main_plot <- ggplot() +
    geom_sf(data = basin_sf, fill = "gray90", color = "gray70", size = 0.5) +
    geom_sf(data = managed_edges_sf, aes(color = production_proportion), size = 1.2) +
    scale_color_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Production\nProportion",
      na.value = "grey60",
      labels = scales::percent_format(accuracy = 1),
      guide = guide_colorbar(
        barwidth = 1, barheight = 15,
        frame.colour = "grey40", ticks.colour = "grey40",
        show.limits = TRUE
      )
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0(subset_label, ": Management Rivers - ", watershed, " Watershed"),
      subtitle = paste("Year", year)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "grey30"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold", color = "grey30"),
      legend.text = element_text(color = "grey30"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  # Apply watershed ordering for bar chart
  mgmt_result_ordered <- apply_watershed_order(mgmt_result, "mgmt_river", reverse_for_plots = TRUE)
  
  # Bar chart
  bar_plot <- ggplot(mgmt_result_ordered, 
                     aes(x = mgmt_river, y = production_proportion)) +
    geom_col(aes(fill = production_proportion), alpha = 0.9) +
    scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), guide = "none") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(title = "Production by Management River", x = "", y = "Production Proportion") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(5, 10, 5, 5, "mm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Plot components
  print(main_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(bar_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  # Add histogram if provided
  if (!is.null(gg_hist)) {
    gg_hist_clean <- gg_hist + 
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    print(gg_hist_clean, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  }
  
  dev.off()
  
  message(paste("Created management map:", basename(output_path)))
}

################################################################################
# AVERAGE PRODUCTION MAPPING FUNCTIONS
################################################################################

#' Create average production map by quartile
create_average_production_map <- function(quartile_data, quartile_label, edges, basin, 
                                          year_range, output_path) {
  
  png(file = output_path, width = 10, height = 8, units = "in", res = 300, bg = "white")
  
  # Join with spatial data
  edges_with_data <- edges %>%
    left_join(quartile_data, by = "mgmt_river") %>%
    filter(!is.na(mgmt_river) & mgmt_river != "") %>%
    mutate(
      production_proportion = ifelse(is.na(production_proportion), 0, production_proportion),
      stream_order = ifelse(is.na(Str_Order), 3, Str_Order),
      line_width = pmax(0.3, pmin(3.0, 0.3 + (stream_order - min(stream_order, na.rm = TRUE)) * 
                                    (3.0 - 0.3) / (max(stream_order, na.rm = TRUE) - min(stream_order, na.rm = TRUE))))
    )
  
  # Ensure consistent CRS
  if (st_crs(basin) != st_crs(edges_with_data)) {
    basin <- st_transform(basin, st_crs(edges_with_data))
  }
  
  # Main map plot
  main_plot <- ggplot() +
    geom_sf(data = basin, fill = "gray95", color = "gray70", 
            linewidth = 0.5, alpha = 0.3) +
    geom_sf(data = edges_with_data, 
            aes(color = production_proportion, linewidth = stream_order), 
            alpha = 0.8) +
    scale_color_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Production\nProportion",
      na.value = "grey60",
      labels = scales::percent_format(accuracy = 1),
      guide = guide_colorbar(
        barwidth = 1, barheight = 15,
        frame.colour = "grey40", ticks.colour = "grey40",
        show.limits = TRUE
      )
    ) +
    scale_linewidth_continuous(
      range = c(0.3, 3.0), 
      name = "Stream\nOrder"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Average ", quartile_label, ": Management Rivers - Kusko Watershed"),
      subtitle = paste("Average across all years (", year_range, ")", sep = "")
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "grey30"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey50"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold", color = "grey30"),
      legend.text = element_text(color = "grey30"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10, "mm")
    )
  
  print(main_plot)
  dev.off()
  
  message(paste("Created average map:", basename(output_path)))
}

#' Create production variability boxplots
create_production_boxplots <- function(mgmt_data, watershed, years, output_dir) {
  
  # Prepare data for boxplots
  quartile_boxplot_data <- mgmt_data %>%
    select(year, mgmt_river, quartile_clean, within_quartile_prop) %>%
    mutate(within_quartile_pct = within_quartile_prop * 100)
  
  # Apply watershed ordering (reverse for coord_flip)
  quartile_boxplot_data <- apply_watershed_order(quartile_boxplot_data, "mgmt_river", reverse_for_plots = TRUE)
  
  year_range <- paste0(min(years), "-", max(years))
  
  # Create individual boxplot for each quartile
  for (q in CONFIG$quartiles) {
    q_data <- quartile_boxplot_data %>% filter(quartile_clean == q)
    
    q_boxplot <- ggplot(q_data, aes(x = mgmt_river, y = within_quartile_pct)) +
      geom_boxplot(
        fill = "grey90", 
        alpha = 0.7, 
        outlier.alpha = 0.8,
        outlier.size = 2.5,
        outlier.shape = 16,
        linewidth = 0.6
      ) +
      scale_y_continuous(
        labels = function(x) paste0(round(x, 1), "%"),
        limits = c(0, max(quartile_boxplot_data$within_quartile_pct, na.rm = TRUE) * 1.05),
        expand = expansion(mult = c(0.02, 0.05))
      ) +
      coord_flip() +
      labs(
        title = paste("Management Unit Production Variability:", q),
        subtitle = paste("Within-", q, " production for each management unit across years (", year_range, ") | Watershed order: upstream → downstream", sep = ""),
        x = "Management Unit (Watershed Position)",
        y = paste(q, "Production (%)"),
        caption = "Boxes show median, quartiles, and whiskers; dots show outliers | Ordered by watershed position (upstream at top)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey50"),
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
    
    ggsave(file.path(output_dir, paste0("management_unit_", q, "_boxplot.png")), 
           q_boxplot, width = 12, height = 10, dpi = 300, bg = "white")
  }
  
  # Create combined faceted boxplot
  combined_boxplot <- ggplot(quartile_boxplot_data, aes(x = mgmt_river, y = within_quartile_pct)) +
    geom_boxplot(
      fill = "grey90", 
      alpha = 0.7, 
      outlier.alpha = 0.8,
      outlier.size = 2,
      outlier.shape = 16,
      linewidth = 0.5
    ) +
    facet_wrap(~quartile_clean, scales = "free_y", ncol = 2) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    coord_flip() +
    labs(
      title = "Management Unit Production Variability Across All Quartiles",
      subtitle = paste("Within-quartile production for each management unit (", year_range, ") | Watershed order: upstream → downstream", sep = ""),
      x = "Management Unit (Watershed Position)",
      y = "Production (%)",
      caption = "Boxes show median, quartiles, and whiskers; dots show outliers | Ordered by watershed position (upstream at top)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey30"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 12),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  ggsave(file.path(output_dir, "management_unit_all_quartiles_boxplot.png"), 
         combined_boxplot, width = 14, height = 12, dpi = 300, bg = "white")
  
  message("Created production variability boxplots")
}

################################################################################
# CUMULATIVE DISTRIBUTION FUNCTIONS
################################################################################

#' Create CPUE histogram for cumulative analysis
create_cpue_histogram <- function(natal_data, year, watershed) {
  doy_to_date <- function(doy, year = 2024) {
    as.Date(doy - 1, origin = paste0(year, "-01-01"))
  }
  
  june_1_doy <- 152
  july_30_doy <- 211
  doy_breaks <- seq(june_1_doy, july_30_doy, by = 10)
  
  ggplot(natal_data, aes(x = DOY, y = dailyCPUEprop)) +
    geom_line(color = "steelblue", linewidth = 1.5, alpha = 0.8) +
    geom_ribbon(aes(ymin = 0, ymax = dailyCPUEprop), 
                fill = "steelblue", alpha = 0.3) +
    scale_x_continuous(
      breaks = doy_breaks,
      labels = function(x) paste0(x, "\n", format(doy_to_date(x), "%b %d")),
      limits = c(june_1_doy, july_30_doy),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(x * 100, 1), "%"),
      limits = c(0, max(natal_data$dailyCPUEprop, na.rm = TRUE) * 1.05)
    ) +
    labs(
      title = paste("Daily CPUE Distribution -", year),
      subtitle = paste(watershed, "Watershed | Fish catch timing throughout the season"),
      x = "Day of Year (Date)",
      y = "Daily CPUE Proportion"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
      axis.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

#' Create cumulative distribution plots
create_cumulative_plots <- function(cumulative_data, output_dir, watershed, interval_days) {
  
  plot_dir <- file.path(output_dir, "Plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get management units and apply watershed ordering
  mgmt_rivers <- sort(unique(cumulative_data$mgmt_river))
  final_order <- WATERSHED_ORDER[WATERSHED_ORDER %in% mgmt_rivers]
  missing_units <- setdiff(mgmt_rivers, WATERSHED_ORDER)
  if (length(missing_units) > 0) {
    final_order <- c(final_order, missing_units)
  }
  
  cumulative_data$mgmt_river <- factor(cumulative_data$mgmt_river, levels = final_order)
  
  # Create colors (red to blue gradient)
  n_units <- length(final_order)
  watershed_colors <- colorRampPalette(c(
    "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
    "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
  ))(n_units)
  names(watershed_colors) <- final_order
  
  # Individual plots for each year
  years <- sort(unique(cumulative_data$year))
  
  for (year in years) {
    year_data <- cumulative_data %>% filter(year == !!year)
    
    # Load natal data for CPUE histogram
    natal_data <- load_natal_data(year, "Kusko")
    cpue_histogram <- create_cpue_histogram(natal_data, year, "Kusko")
    
    # Create timing plot
    p_year <- ggplot(year_data, aes(x = doy, y = cumulative_percent, color = mgmt_river)) +
      geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.7) +
      geom_hline(yintercept = c(10, 90), color = "gray80", linetype = "dotted", alpha = 0.5) +
      geom_line(linewidth = 2.5, alpha = 0.9) +
      scale_color_manual(values = watershed_colors, 
                         breaks = final_order,
                         guide = guide_legend(override.aes = list(linewidth = 3))) +
      scale_y_continuous(
        labels = function(x) paste0(round(x, 1), "%"),
        limits = c(0, 100),
        breaks = c(0, 10, 25, 50, 75, 90, 100),
        minor_breaks = seq(0, 100, 5)
      ) +
      scale_x_continuous(
        breaks = sort(unique(year_data$doy)),
        labels = function(x) {
          date_str <- format(as.Date(x - 1, origin = paste0(year, "-01-01")), "%b %d")
          paste0("DOY ", x, "\n", date_str)
        },
        limits = c(min(year_data$doy), max(year_data$doy)),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      labs(
        title = paste("Salmon Run Timing Progress -", year),
        subtitle = paste("Kusko Watershed | Management unit progression |", interval_days, "day intervals"),
        x = paste("Day of Year (Date) -", interval_days, "Day Intervals"),
        y = "Cumulative Percent of Management Unit's Total Production",
        color = "Management Unit\n(Watershed Position)",
        caption = paste("Red=upstream, Blue=downstream |", length(unique(year_data$doy)), "data points")
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "gray20"),
        plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 10, hjust = 0.5, color = "gray50"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        axis.line = element_line(color = "gray20", linewidth = 1.2),
        axis.text = element_text(size = 10, color = "gray20"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, color = "gray20"),
        axis.title = element_text(face = "bold", size = 12, color = "gray20"),
        axis.ticks = element_line(color = "gray20", linewidth = 0.8),
        panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
        panel.border = element_rect(color = "gray20", fill = NA, linewidth = 1.0),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    # [trimmed: not a paper figure] per-year run-timing progress plots not exported.
    # combined_filename <- paste0("run_timing_progress_", year, "_", watershed, ".png")
    # png(file.path(plot_dir, combined_filename), width = 18, height = 12, units = "in", res = 300, bg = "white")
    # grid.arrange(p_year, cpue_histogram, nrow = 2, heights = c(3, 1))
    # dev.off()
    #
    # message(paste("Created timing plot for", year))
  }
  
  # Create overview faceted plot
  doy_range_overview <- range(cumulative_data$doy)
  
  p1_faceted <- ggplot(cumulative_data, aes(x = doy, y = cumulative_percent, color = mgmt_river)) +
    geom_line(linewidth = 2.0, alpha = 0.9) +
    facet_wrap(~year, ncol = 1, scales = "free_x") +
    scale_color_manual(values = watershed_colors,
                       breaks = final_order,
                       guide = guide_legend(override.aes = list(linewidth = 3),
                                            ncol = 1,
                                            keywidth = 1.5,
                                            keyheight = 1.0)) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, 100),
      breaks = seq(0, 100, 25)
    ) +
    scale_x_continuous(
      limits = doy_range_overview,
      breaks = seq(from = min(cumulative_data$doy), to = max(cumulative_data$doy), by = 9),
      labels = function(x) {
        date_obj <- as.Date(x - 1, origin = "2024-01-01")
        format(date_obj, "%b %d")
      },
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = "Cumulative Distribution of Returning Chinook Stocks Over Time",
      x = "Date",
      y = "Cumulative Percent of Unit's Total Production",
      color = "Management Unit"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 11),
      legend.margin = margin(l = 20),
      legend.box.spacing = unit(0.5, "cm"),
      axis.line = element_line(color = "gray20", linewidth = 1.0),
      axis.text = element_text(color = "gray20", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "gray20"),
      axis.title = element_text(face = "bold", color = "gray20", size = 13),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.4),
      panel.spacing = unit(1.5, "lines"),
      plot.margin = margin(15, 15, 15, 15, "mm")
    )
  
  # [trimmed: not the published Fig 4] management-unit cumulative overview is not
  # exported. The published Fig 4 is the regional-group MULTIPANEL from R/06; the
  # cumulative-distribution CSV intermediate is still written by R/03.
  # ggsave(file.path(plot_dir, paste0("cumulative_overview_", watershed, "_PAPER.png")),
  #        p1_faceted, width = 16, height = 20, dpi = 300, bg = "white")
  
  message("Created cumulative distribution plots")
}

################################################################################
# CLOSURE PROTECTION FUNCTIONS
################################################################################

#' Create closure protection boxplot
create_closure_boxplot <- function(closure_data, output_dir) {
  
  year_range <- paste0(min(closure_data$year), "-", max(closure_data$year))
  
  # Get management units and apply watershed ordering (same as cumulative)
  mgmt_rivers <- sort(unique(closure_data$mgmt_river))
  final_order <- WATERSHED_ORDER[WATERSHED_ORDER %in% mgmt_rivers]
  missing_units <- setdiff(mgmt_rivers, WATERSHED_ORDER)
  if (length(missing_units) > 0) {
    final_order <- c(final_order, missing_units)
  }
  
  # Apply ordering (reverse for coord_flip so upstream appears at top)
  closure_data$mgmt_river <- factor(closure_data$mgmt_river, levels = rev(final_order))
  
  # Create colors EXACTLY matching cumulative distribution plots
  n_units <- length(final_order)
  watershed_colors <- colorRampPalette(c(
    "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
    "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
  ))(n_units)
  names(watershed_colors) <- final_order
  
  # Reverse the color mapping to match the reversed factor levels
  watershed_colors_reversed <- watershed_colors[rev(names(watershed_colors))]
  names(watershed_colors_reversed) <- rev(final_order)
  
  # CONTEMPORARY DESIGN: Range plot with distributions
  boxplot <- ggplot(closure_data, aes(y = mgmt_river, x = protection_percentage, color = mgmt_river, fill = mgmt_river)) +
    
    # Modern approach: violin plot for distribution shape (subtle)
    geom_violin(
      alpha = 0.15,
      scale = "width",
      width = 0.6,
      color = NA  # No violin borders
    ) +
    
    # Add subtle range lines (min to max)
    stat_summary(
      fun.min = min, 
      fun.max = max,
      geom = "linerange",
      linewidth = 2,
      alpha = 0.6
    ) +
    
    # Prominent MEAN points (changed from median)
    stat_summary(
      fun = mean,  # Changed to mean
      geom = "point",
      size = 4,
      color = "white",
      stroke = 2,
      shape = 21  # Circle with border
    ) +
    
    
    
    # Modern color scheme
    scale_color_manual(values = watershed_colors_reversed, guide = "none") +
    scale_fill_manual(values = watershed_colors_reversed, guide = "none") +
    
    # Clean modern axes
    scale_x_continuous(
      name = "Protection effectiveness (%)",
      breaks = function(x) pretty(x, n = 6),
      expand = expansion(mult = c(0.02, 0.05)),
      labels = function(x) paste0(round(x, 1), "%")
    ) +
    
    scale_y_discrete(
      name = NULL,  # Remove y-axis title for cleaner look
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    # Contemporary minimal title
    labs(
      title = "Front-End Closure Protection Effectiveness",
      subtitle = paste("Distribution across", year_range, "• Red = Upstream • Blue = Downstream"),
      caption = "White dots = mean • Colored areas show data distribution • Range lines show min-max"
    ) +
    
    # Ultra-modern theme
    theme_void() +
    theme(
      # Contemporary typography
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0,
        color = "gray15",
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = 14,
        hjust = 0,
        color = "gray50",
        margin = margin(b = 25)
      ),
      plot.caption = element_text(
        size = 11,
        hjust = 0,
        color = "gray60",
        margin = margin(t = 20)
      ),
      
      # Minimal axis text
      axis.text.y = element_text(
        size = 12,
        color = "gray30",
        hjust = 1,
        margin = margin(r = 15)
      ),
      axis.text.x = element_text(
        size = 11,
        color = "gray30",
        margin = margin(t = 10)
      ),
      axis.title.x = element_text(
        size = 13,
        color = "gray20",
        face = "bold",
        margin = margin(t = 15)
      ),
      
      # Extremely minimal grid
      panel.grid.major.x = element_line(
        color = "gray94",
        linewidth = 0.3
      ),
      panel.grid.minor.x = element_line(
        color = "gray97",
        linewidth = 0.2
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      # Clean spacing
      plot.margin = margin(30, 30, 30, 30, "mm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(output_dir, "front_end_closure_protection_modern.png"), 
         boxplot, width = 14, height = 10, dpi = 300, bg = "white")
  
  message("Created contemporary closure protection plot with modern design")
  message("Features: violin distributions + range lines + median points + individual data")
  
  # Print color mapping for verification
  message("Color mapping (upstream to downstream):")
  for (i in 1:length(final_order)) {
    message(paste("  ", final_order[i], ":", watershed_colors[i]))
  }
}

### Traditional boxplot 
create_closure_boxplot_traditional <- function(closure_data, output_dir) {
  
  year_range <- paste0(min(closure_data$year), "-", max(closure_data$year))
  
  # Get management units and apply watershed ordering (same as cumulative)
  mgmt_rivers <- sort(unique(closure_data$mgmt_river))
  final_order <- WATERSHED_ORDER[WATERSHED_ORDER %in% mgmt_rivers]
  missing_units <- setdiff(mgmt_rivers, WATERSHED_ORDER)
  if (length(missing_units) > 0) {
    final_order <- c(final_order, missing_units)
  }
  
  # Apply ordering (reverse for coord_flip so upstream appears at top)
  closure_data$mgmt_river <- factor(closure_data$mgmt_river, levels = rev(final_order))
  
  # Create colors EXACTLY matching cumulative distribution plots
  n_units <- length(final_order)
  watershed_colors <- colorRampPalette(c(
    "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
    "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
  ))(n_units)
  names(watershed_colors) <- final_order
  
  # Reverse the color mapping to match the reversed factor levels
  watershed_colors_reversed <- watershed_colors[rev(names(watershed_colors))]
  names(watershed_colors_reversed) <- rev(final_order)
  
  # MODERN TRADITIONAL BOXPLOT - KEEPING ORIGINAL FILL COLORS
  boxplot <- ggplot(closure_data, aes(x = mgmt_river, y = protection_percentage, fill = mgmt_river)) +
    geom_boxplot(
      alpha = 0.8,  # Keep original alpha
      outlier.size = 2.8,  # Keep original size
      outlier.shape = 16,
      outlier.alpha = 0.9,
      outlier.colour = "gray20",  # Keep original outlier color
      colour = "gray30",  # Keep original box border color
      linewidth = 0.5,  # Keep original line width
      width = 0.75  # Keep original width
    ) +
    
    # Use the exact same colors as cumulative distribution plots - UNCHANGED
    scale_fill_manual(values = watershed_colors_reversed, guide = "none") +
    coord_flip() +
    
    # Modern axis styling 
    scale_y_continuous(
      name = "Protection Effectiveness (%)",
      breaks = function(x) pretty(x, n = 6),
      expand = expansion(mult = c(0.02, 0.05)),
      labels = function(x) paste0(round(x, 1), "%")
    ) +
    
    scale_x_discrete(name = "Management Unit") +
    
    # Clean, modern title - NO SUBTITLES
    labs(
      title = "Front-End Closure Protection by Management Unit"
    ) +
    
    # Modern minimal theme
    theme_minimal(base_size = 14) +
    theme(
      # Modern typography
      plot.title = element_text(
        family = "Arial",
        face = "bold", 
        size = 20, 
        hjust = 0.5, 
        color = "#2c3e50",
        margin = margin(b = 25)
      ),
      
      # Modern axis styling - NO SUBTITLE TEXT
      axis.text.y = element_text(
        family = "Arial",
        size = 12, 
        color = "#2c3e50"
      ),
      axis.text.x = element_text(
        family = "Arial",
        size = 11, 
        color = "#2c3e50"
      ),
      axis.title = element_text(
        family = "Arial",
        face = "bold", 
        size = 14, 
        color = "#34495e"
      ),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      
      # Ultra-minimal grid
      panel.grid.major.x = element_line(color = "#ecf0f1", linewidth = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      # Modern clean spacing
      plot.margin = margin(25, 25, 25, 25, "mm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Remove traditional border for cleaner look
      panel.border = element_blank(),
      
      # Add subtle axis lines
      axis.line.x = element_line(color = "#bdc3c7", linewidth = 0.5),
      axis.ticks = element_line(color = "#bdc3c7", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm")
    )
  
  ggsave(file.path(output_dir, "front_end_closure_protection_boxplot_modern.png"), 
         boxplot, width = 12, height = 10, dpi = 300, bg = "white")
  
  message("Created modern traditional closure protection boxplot")
}

# Master function that creates BOTH plots
create_closure_plots <- function(closure_data, output_dir) {
  
  # Create the modern violin/range plot
  create_closure_boxplot(closure_data, output_dir)
  
  # Create the traditional boxplot
  create_closure_boxplot_traditional(closure_data, output_dir)
  
  message("Created both modern and traditional closure protection plots:")
  message("  - front_end_closure_protection_modern.png (violin + mean)")
  message("  - front_end_closure_protection_boxplot_traditional.png (classic boxplot)")
}


cat("✓ All visualization functions loaded successfully\n")