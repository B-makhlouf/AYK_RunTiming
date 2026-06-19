################################################################################
# ISOTOPE CLUSTERING ANALYSIS - ELBOW METHOD
# Using ALL isotope values from ALL management units
# EXCLUDING: Johnson and Lower Kusko
# COMBINING: South Fork Kusko INTO Upper Kuskokwim
################################################################################
# This script runs automatically when sourced
################################################################################

run_isotope_clustering <- function() {
  
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(cluster)
  library(factoextra)
  library(gridExtra)
  library(viridis)
  
  # ============================================================================
  # CONFIGURATION
  # ============================================================================
  
  # File paths come from config.R (single source of truth)
  if (!exists("PROJECT_ROOT")) source("config.R")
  EDGES_PATH <- KUSKO_EDGES
  BASIN_PATH <- KUSKO_BASIN

  # Output directory
  OUTPUT_DIR <- CLUSTER_GROUPS_DIR
  
  # Management unit ordering (upstream to downstream)
  # EXCLUDING: Johnson and Lower Kusko
  # COMBINING: S. Fork Kusko merged into Upper Kusko Main
  MGMT_ORDER <- c(
    "N. Fork Kusko", 
    "E. Fork Kuskokwim",
    "Upper Kusko Main",  # Now includes S. Fork Kusko data
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
    "Kwethluk"
    # S. Fork Kusko - COMBINED with Upper Kusko Main
    # Johnson - EXCLUDED
    # Lower Kusko - EXCLUDED
  )
  
  # ============================================================================
  # LOAD AND PREPARE DATA
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("LOADING ISOTOPE DATA\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Load edges shapefile
  edges <- st_read(EDGES_PATH, quiet = TRUE)
  
  # Filter to edges with management unit assignments and isotope data
  # Apply same filters as in your analysis: stream order, priors, etc.
  managed_edges <- edges %>%
    filter(
      !is.na(mgmt_river), 
      mgmt_river != "", 
      !is.na(iso_pred),
      Str_Order >= 3,  # Minimum stream order from your setup
      UniPh2oNoE > 0   # Applying Kusko prior
    )
  
  # Additional filter: only spawning habitat (matching your PresencePrior)
  managed_edges <- managed_edges %>%
    filter(!(Str_Order %in% c(6, 7, 8) & SPAWNING_C == 0))
  
  # EXCLUDE Johnson and Lower Kusko from analysis
  managed_edges <- managed_edges %>%
    filter(!mgmt_river %in% c("Johnson", "Lower Kusko"))
  
  # COMBINE S. Fork Kusko into Upper Kusko Main
  managed_edges <- managed_edges %>%
    mutate(mgmt_river = case_when(
      mgmt_river == "S. Fork Kusko" ~ "Upper Kusko Main",
      TRUE ~ mgmt_river
    ))
  
  cat(sprintf("Loaded %d stream edges with isotope data\n", nrow(managed_edges)))
  cat(sprintf("Covering %d management units\n", length(unique(managed_edges$mgmt_river))))
  cat("\nCOMBINATIONS/EXCLUSIONS APPLIED:\n")
  cat("  - S. Fork Kusko: COMBINED into Upper Kusko Main\n")
  cat("  - Johnson: EXCLUDED from analysis\n")
  cat("  - Lower Kusko: EXCLUDED from analysis\n\n")
  
  # ============================================================================
  # EXTRACT ALL ISOTOPE VALUES (NOT JUST MEANS!)
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("EXTRACTING ALL ISOTOPE VALUES FROM ALL MANAGEMENT UNITS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Extract ALL isotope values with their management units
  # THIS IS THE KEY: We use every single iso_pred value, not aggregated means
  isotope_data <- managed_edges %>%
    st_drop_geometry() %>%
    select(mgmt_river, iso_pred) %>%
    filter(!is.na(iso_pred))
  
  cat(sprintf("Total isotope values extracted: %d\n", nrow(isotope_data)))
  cat("These are ALL individual isotope predictions from all stream edges,\n")
  cat("NOT averaged or aggregated by management unit.\n")
  cat("S. Fork Kusko values are now included in Upper Kusko Main.\n\n")
  
  # Summary statistics by management unit (for reporting)
  mgmt_summary <- isotope_data %>%
    group_by(mgmt_river) %>%
    summarise(
      mean_iso = mean(iso_pred, na.rm = TRUE),
      sd_iso = sd(iso_pred, na.rm = TRUE),
      min_iso = min(iso_pred, na.rm = TRUE),
      max_iso = max(iso_pred, na.rm = TRUE),
      n_edges = n(),
      .groups = "drop"
    ) %>%
    arrange(match(mgmt_river, MGMT_ORDER))
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("ISOTOPE SUMMARY BY MANAGEMENT UNIT\n")
  cat("(Upper Kusko Main now includes S. Fork Kusko data)\n")
  cat(rep("=", 80), "\n", sep = "")
  print(as.data.frame(mgmt_summary), row.names = FALSE)
  
  # ============================================================================
  # PREPARE DATA FOR CLUSTERING - USING ALL VALUES
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("PREPARING DATA FOR CLUSTERING\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Use ALL isotope values (not just means) for clustering
  # This captures the full distribution within each management unit
  cluster_data <- isotope_data %>%
    arrange(match(mgmt_river, MGMT_ORDER))
  
  # Create matrix for clustering - ALL isotope values
  X <- as.matrix(cluster_data$iso_pred)
  
  # Create labels showing which management unit each value belongs to
  mgmt_labels <- cluster_data$mgmt_river
  
  # Standardize the data
  X_scaled <- scale(X)
  
  cat(sprintf("✓ Using ALL %d isotope values for clustering\n", nrow(X)))
  cat(sprintf("✓ From %d management units\n", length(unique(mgmt_labels))))
  cat(sprintf("✓ Isotope range: %.4f to %.4f\n", min(X), max(X)))
  cat(sprintf("✓ Mean: %.4f, SD: %.4f\n", mean(X), sd(X)))
  cat("\nValues per management unit:\n")
  cat("(Note: Upper Kusko Main counts include former S. Fork Kusko)\n")
  print(table(mgmt_labels))
  
  # ============================================================================
  # ELBOW METHOD ANALYSIS
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("PERFORMING ELBOW METHOD ANALYSIS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Test different numbers of clusters
  K_range <- 2:10
  wss <- numeric(length(K_range))  # Within-cluster sum of squares
  silhouette_avg <- numeric(length(K_range))
  between_ss <- numeric(length(K_range))
  
  for (i in seq_along(K_range)) {
    k <- K_range[i]
    set.seed(42)
    kmeans_result <- kmeans(X_scaled, centers = k, nstart = 25)
    wss[i] <- kmeans_result$tot.withinss
    between_ss[i] <- kmeans_result$betweenss
    
    # Calculate silhouette score
    if (k > 1 && k < nrow(X_scaled)) {
      sil <- silhouette(kmeans_result$cluster, dist(X_scaled))
      silhouette_avg[i] <- mean(sil[, 3])
    } else {
      silhouette_avg[i] <- NA
    }
    
    cat(sprintf("k = %2d | WSS = %.4f | Between SS = %.4f | Silhouette = %.4f\n", 
                k, wss[i], between_ss[i], silhouette_avg[i]))
  }
  
  # Calculate rate of change in WSS
  wss_change <- c(NA, -diff(wss))
  wss_pct_change <- c(NA, (-diff(wss) / wss[-length(wss)]) * 100)
  
  # ============================================================================
  # VISUALIZATION
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("CREATING VISUALIZATIONS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Create output directory
  output_dir <- OUTPUT_DIR
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- PLOT 1: Elbow Plot ---
  p1 <- ggplot(data.frame(k = K_range, wss = wss), aes(x = k, y = wss)) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_point(color = "steelblue", size = 3) +
    geom_vline(xintercept = which.max(diff(diff(wss))) + 2, 
               linetype = "dashed", color = "red", alpha = 0.5) +
    labs(
      title = "Elbow Method",
      subtitle = "Within-cluster Sum of Squares",
      x = "Number of Clusters (k)",
      y = "Total WSS"
    ) +
    scale_x_continuous(breaks = K_range) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey50"),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # --- PLOT 2: Silhouette Score ---
  p2 <- ggplot(data.frame(k = K_range, sil = silhouette_avg), 
               aes(x = k, y = sil)) +
    geom_line(color = "darkgreen", linewidth = 1.2) +
    geom_point(color = "darkgreen", size = 3) +
    geom_vline(xintercept = K_range[which.max(silhouette_avg)], 
               linetype = "dashed", color = "red", alpha = 0.5) +
    labs(
      title = "Silhouette Analysis",
      subtitle = "Average Silhouette Width",
      x = "Number of Clusters (k)",
      y = "Avg Silhouette Score"
    ) +
    scale_x_continuous(breaks = K_range) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey50"),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # --- PLOT 3: Variance Explained ---
  total_ss <- sum((X_scaled - mean(X_scaled))^2)
  variance_explained <- (between_ss / total_ss) * 100
  
  p3 <- ggplot(data.frame(k = K_range, var_exp = variance_explained), 
               aes(x = k, y = var_exp)) +
    geom_line(color = "purple", linewidth = 1.2) +
    geom_point(color = "purple", size = 3) +
    labs(
      title = "Variance Explained",
      subtitle = "Between-cluster variance",
      x = "Number of Clusters (k)",
      y = "Variance Explained (%)"
    ) +
    scale_x_continuous(breaks = K_range) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey50"),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # --- PLOT 4: Rate of WSS Decrease ---
  p4 <- ggplot(data.frame(k = K_range, decrease = wss_pct_change), 
               aes(x = k, y = decrease)) +
    geom_line(color = "orange", linewidth = 1.2, na.rm = TRUE) +
    geom_point(color = "orange", size = 3, na.rm = TRUE) +
    labs(
      title = "Rate of Improvement",
      subtitle = "% Decrease in WSS",
      x = "Number of Clusters (k)",
      y = "% Decrease in WSS"
    ) +
    scale_x_continuous(breaks = K_range) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey50"),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                                top = grid::textGrob(
                                  "Elbow Method Analysis - All Isotope Values\nUpper Kusko includes S. Fork Kusko",
                                  gp = grid::gpar(fontsize = 14, fontface = "bold")
                                ))
  
  # Save combined plot
  ggsave(
    filename = file.path(output_dir, "elbow_method_analysis_combined_upper_kusko.png"),
    plot = combined_plot,
    width = 12,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  cat(sprintf("✓ Saved 4-panel elbow analysis plot to: %s\n", 
              file.path(output_dir, "elbow_method_analysis_combined_upper_kusko.png")))
  
  # ============================================================================
  # DETAILED CLUSTERING FOR OPTIMAL K VALUES
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("PERFORMING DETAILED CLUSTERING FOR KEY K VALUES\n")
  cat(rep("=", 80), "\n", sep = "")
  
  optimal_k_values <- c(3, 4, 5, 6)  # Test multiple k values
  cluster_results <- list()
  
  for (k in optimal_k_values) {
    cat(sprintf("\n--- Clustering with k = %d ---\n", k))
    
    set.seed(42)
    kmeans_result <- kmeans(X_scaled, centers = k, nstart = 25)
    
    # Create data frame with ALL values and their cluster assignments
    clustered_data <- data.frame(
      mgmt_river = mgmt_labels,
      iso_pred = X[, 1],
      cluster = kmeans_result$cluster
    )
    
    # Calculate cluster statistics
    cluster_stats <- clustered_data %>%
      group_by(cluster) %>%
      summarise(
        mean_iso = mean(iso_pred),
        sd_iso = sd(iso_pred),
        min_iso = min(iso_pred),
        max_iso = max(iso_pred),
        n_values = n(),
        .groups = "drop"
      ) %>%
      arrange(mean_iso)
    
    # Get management unit composition of each cluster
    cluster_composition <- clustered_data %>%
      group_by(cluster, mgmt_river) %>%
      summarise(n_values = n(), .groups = "drop") %>%
      group_by(cluster) %>%
      mutate(
        pct_of_cluster = (n_values / sum(n_values)) * 100
      ) %>%
      arrange(cluster, desc(n_values))
    
    # For summarizing by management unit (what cluster is each unit primarily in?)
    cluster_assignment <- clustered_data %>%
      group_by(mgmt_river) %>%
      summarise(
        cluster = as.integer(names(sort(table(cluster), decreasing = TRUE)[1])),
        pct_in_cluster = (max(table(cluster)) / n()) * 100,
        mean_iso = mean(iso_pred),
        n_values = n(),
        .groups = "drop"
      ) %>%
      arrange(match(mgmt_river, MGMT_ORDER))
    
    cat("\nCluster Statistics (ordered by mean isotope):\n")
    print(as.data.frame(cluster_stats), row.names = FALSE)
    
    cat("\nManagement Unit Cluster Assignments:\n")
    cat("(Upper Kusko Main includes S. Fork Kusko data)\n")
    print(as.data.frame(cluster_assignment), row.names = FALSE)
    
    # Store results
    cluster_results[[as.character(k)]] <- cluster_assignment
    
    # --- Visualization: Distribution by cluster ---
    # Order clusters by mean isotope value
    clustered_data$cluster_ordered <- factor(
      clustered_data$cluster,
      levels = cluster_stats$cluster,
      labels = paste0("Cluster ", cluster_stats$cluster, "\n(μ=", 
                      round(cluster_stats$mean_iso, 3), ")")
    )
    
    dist_plot <- ggplot(clustered_data, 
                        aes(x = iso_pred, fill = cluster_ordered)) +
      geom_histogram(bins = 50, alpha = 0.7, color = "white", linewidth = 0.3) +
      scale_fill_viridis_d(option = "D") +
      labs(
        title = sprintf("Distribution of Isotope Values by Cluster (k=%d)", k),
        subtitle = sprintf("Using ALL %d isotope values - Upper Kusko includes S. Fork Kusko",
                           nrow(isotope_data)),
        x = "Isotope Value (δ²H)",
        y = "Count",
        fill = "Cluster"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey50"),
        legend.position = "right"
      )
    
    ggsave(
      filename = file.path(output_dir, 
                           sprintf("clusters_k%d_combined_upper_kusko.png", k)),
      plot = dist_plot,
      width = 12,
      height = 6,
      dpi = 300,
      bg = "white"
    )
    
    cat(sprintf("✓ Saved distribution plot to: %s\n", 
                file.path(output_dir, sprintf("clusters_k%d_combined_upper_kusko.png", k))))
    
    # --- Spatial Visualization: Basin Map ---
    if (file.exists(BASIN_PATH)) {
      cat(sprintf("Creating basin map for k=%d...\n", k))
      
      # Load basin
      basin <- st_read(BASIN_PATH, quiet = TRUE)
      
      # Join cluster assignments to edges
      # Since we combined S. Fork into Upper Kusko, the join will work automatically
      edges_with_clusters <- edges %>%
        mutate(mgmt_river_adjusted = case_when(
          mgmt_river == "S. Fork Kusko" ~ "Upper Kusko Main",
          TRUE ~ mgmt_river
        )) %>%
        left_join(cluster_assignment %>% select(mgmt_river, cluster), 
                  by = c("mgmt_river_adjusted" = "mgmt_river")) %>%
        mutate(
          cluster = as.character(cluster),
          stream_order = pmin(Str_Order, 8)  # Cap at 8 for visualization
        )
      
      # Split into clustered and non-clustered edges
      clustered_edges <- edges_with_clusters %>% filter(!is.na(cluster))
      non_clustered_edges <- edges_with_clusters %>% filter(is.na(cluster) | cluster == "")
      
      cat(sprintf("Plotting %d clustered edges and %d non-clustered edges\n", 
                  nrow(clustered_edges), nrow(non_clustered_edges)))
      
      # Get cluster order by mean isotope value for legend
      cluster_order_for_legend <- cluster_assignment %>%
        group_by(cluster) %>%
        summarise(mean_iso = mean(mean_iso), .groups = "drop") %>%
        arrange(mean_iso)
      
      # Create factor with ordered levels for legend
      if (k == 3) {
        legend_labels <- paste("Cluster", 1:3)
      } else if (k == 4) {
        legend_labels <- paste("Cluster", 1:4)
      } else if (k == 5) {
        legend_labels <- paste("Cluster", 1:5)
      } else if (k == 6) {
        legend_labels <- paste("Cluster", 1:6)
      } else {
        legend_labels <- paste("Cluster", 1:k)
      }
      
      edges_with_clusters_labeled <- clustered_edges %>%
        mutate(cluster_factor = factor(cluster, 
                                       levels = cluster_order_for_legend$cluster,
                                       labels = legend_labels))
      
      # Create the map
      map_plot <- ggplot() +
        geom_sf(data = basin, fill = "gray95", color = "gray70", 
                linewidth = 0.5, alpha = 0.3) +
        geom_sf(data = non_clustered_edges,
                aes(linewidth = stream_order),
                color = "gray75",
                alpha = 0.4) +
        geom_sf(data = edges_with_clusters_labeled, 
                aes(color = cluster_factor, linewidth = stream_order), 
                alpha = 0.8) +
        scale_color_viridis_d(
          option = "D",
          name = "Cluster",
          na.value = "grey60"
        ) +
        scale_linewidth_continuous(
          range = c(0.3, 3.0), 
          name = "Stream\nOrder"
        ) +
        coord_sf(datum = NA) +
        labs(
          title = sprintf("Management Units Clustered by Isotope Values (k=%d)", k),
          subtitle = sprintf("Using ALL %d isotope values - Upper Kusko includes S. Fork Kusko",
                             nrow(isotope_data))
        ) +
        theme_void() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "grey30"),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
          legend.position = "right",
          legend.title = element_text(size = 11, face = "bold", color = "grey30"),
          legend.text = element_text(size = 10, color = "grey30"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(10, 10, 10, 10, "mm")
        )
      
      ggsave(
        filename = file.path(output_dir, sprintf("basin_map_k%d_combined_upper_kusko.png", k)),
        plot = map_plot,
        width = 12,
        height = 10,
        dpi = 300,
        bg = "white"
      )
      
      cat(sprintf("✓ Saved basin map to: %s\n", 
                  file.path(output_dir, sprintf("basin_map_k%d_combined_upper_kusko.png", k))))
      
    } else {
      cat(sprintf("⚠ Basin shapefile not found. Skipping basin map creation\n"))
    }
    
    # --- Violin Plot for k=6 ONLY ---
    if (k == 6) {
      cat(sprintf("\n--- Creating Publication-Quality Violin Plot for k=6 ---\n"))
      
      # Prepare data for violin plot
      clustered_data_violin <- clustered_data %>%
        mutate(cluster_factor = factor(cluster))
      
      # Order clusters by mean isotope value and create numeric labels
      cluster_order_violin <- cluster_stats %>%
        arrange(mean_iso) %>%
        mutate(
          cluster_label = paste0("Cluster ", 1:6)
        )
      
      # Join labels to data
      clustered_data_violin <- clustered_data_violin %>%
        left_join(cluster_order_violin %>% select(cluster, cluster_label), by = "cluster") %>%
        mutate(cluster_label = factor(cluster_label, 
                                      levels = paste0("Cluster ", 1:6)))
      
      # Define modern colorblind-friendly palette (viridis)
      violin_colors <- viridis(6, option = "D", begin = 0.15, end = 0.95)
      
      # Create publication-quality violin plot
      p_violin <- ggplot(clustered_data_violin, 
                         aes(x = cluster_label, y = iso_pred, fill = cluster_label)) +
        
        # Violin plot with semi-transparency
        geom_violin(alpha = 0.85, scale = "width", trim = FALSE, 
                    color = "white", linewidth = 1.2) +
        
        # Boxplot overlay showing quartiles
        geom_boxplot(width = 0.12, alpha = 0.95, outlier.shape = NA,
                     color = "gray10", fill = "white", linewidth = 0.6) +
        
        # Median points for emphasis
        stat_summary(fun = median, geom = "point", 
                     size = 3, color = "gray10", shape = 21, 
                     fill = "white", stroke = 1.5) +
        
        # Color scale
        scale_fill_manual(values = violin_colors) +
        
        # Labels with proper strontium notation
        labs(
          title = "Isotope Distribution by Cluster (k=6)",
          subtitle = "Kuskokwim Basin Management Units • Upper Kuskokwim includes S. Fork Kusko",
          x = NULL,
          y = expression(paste(""^"87", "Sr/", ""^"86", "Sr"))
        ) +
        
        # Modern, clean theme
        theme_minimal(base_size = 13, base_family = "sans") +
        theme(
          # Remove legend (redundant with x-axis)
          legend.position = "none",
          
          # Plot titles - left aligned for modern look
          plot.title = element_text(size = 18, face = "bold", 
                                    hjust = 0, color = "#1a1a1a",
                                    margin = margin(b = 8)),
          plot.subtitle = element_text(size = 12, color = "#666666",
                                       hjust = 0, margin = margin(b = 20)),
          
          # Axes
          axis.title.y = element_text(size = 13, face = "bold", 
                                      color = "#333333", margin = margin(r = 12)),
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, 
                                     vjust = 1, color = "#333333", face = "bold"),
          axis.text.y = element_text(size = 11, color = "#4d4d4d"),
          axis.line.x = element_line(color = "#333333", linewidth = 0.6),
          axis.line.y = element_line(color = "#333333", linewidth = 0.6),
          axis.ticks = element_line(color = "#333333", linewidth = 0.6),
          axis.ticks.length = unit(0.15, "cm"),
          
          # Grid - subtle horizontal lines only
          panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.4),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          
          # Background
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          
          # Margins
          plot.margin = margin(25, 25, 20, 25)
        )
      
      # Save clean version
      ggsave(
        filename = file.path(output_dir, "isotope_violin_plot_k6_publication.png"),
        plot = p_violin,
        width = 12,
        height = 7,
        dpi = 600,
        bg = "white"
      )
      
      cat(sprintf("✓ Saved violin plot to: %s\n", 
                  file.path(output_dir, "isotope_violin_plot_k6_publication.png")))
      
      # Create annotated version with sample sizes
      y_min <- min(clustered_data_violin$iso_pred)
      y_range <- max(clustered_data_violin$iso_pred) - y_min
      y_pos <- y_min - (y_range * 0.08)
      
      p_violin_annotated <- p_violin +
        geom_text(data = cluster_order_violin,
                  aes(x = cluster_label, y = y_pos, label = paste0("n = ", formatC(n_values, format="d", big.mark=","))),
                  inherit.aes = FALSE,
                  size = 3.8, color = "#666666", fontface = "italic") +
        scale_y_continuous(expand = expansion(mult = c(0.12, 0.05)))
      
      # Save annotated version
      ggsave(
        filename = file.path(output_dir, "isotope_violin_plot_k6_annotated.png"),
        plot = p_violin_annotated,
        width = 12,
        height = 7,
        dpi = 600,
        bg = "white"
      )
      
      cat(sprintf("✓ Saved annotated violin plot to: %s\n", 
                  file.path(output_dir, "isotope_violin_plot_k6_annotated.png")))
      
      # Create horizontal version (alternative layout)
      p_violin_horizontal <- ggplot(clustered_data_violin, 
                                    aes(x = iso_pred, y = cluster_label, fill = cluster_label)) +
        
        # Horizontal violin plot
        geom_violin(alpha = 0.85, scale = "width", trim = FALSE, 
                    color = "white", linewidth = 1.2) +
        
        # Boxplot overlay
        geom_boxplot(width = 0.12, alpha = 0.95, outlier.shape = NA,
                     color = "gray10", fill = "white", linewidth = 0.6) +
        
        # Median points
        stat_summary(fun = median, geom = "point", 
                     size = 3, color = "gray10", shape = 21, 
                     fill = "white", stroke = 1.5) +
        
        # Color scale - reversed for horizontal
        scale_fill_manual(values = rev(violin_colors)) +
        
        # Reverse y-axis so Extreme Upstream is on top
        scale_y_discrete(limits = rev(levels(clustered_data_violin$cluster_label))) +
        
        # Labels
        labs(
          title = "Isotope Distribution by Cluster (k=6)",
          subtitle = "Kuskokwim Basin Management Units • Upper Kuskokwim includes S. Fork Kusko",
          y = NULL,
          x = expression(paste(""^"87", "Sr/", ""^"86", "Sr"))
        ) +
        
        # Theme
        theme_minimal(base_size = 13, base_family = "sans") +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 18, face = "bold", 
                                    hjust = 0, color = "#1a1a1a",
                                    margin = margin(b = 8)),
          plot.subtitle = element_text(size = 12, color = "#666666",
                                       hjust = 0, margin = margin(b = 20)),
          axis.title.x = element_text(size = 13, face = "bold", 
                                      color = "#333333", margin = margin(t = 12)),
          axis.text.y = element_text(size = 11, color = "#333333", face = "bold",
                                     hjust = 1),
          axis.text.x = element_text(size = 11, color = "#4d4d4d"),
          axis.line = element_line(color = "#333333", linewidth = 0.6),
          axis.ticks = element_line(color = "#333333", linewidth = 0.6),
          axis.ticks.length = unit(0.15, "cm"),
          panel.grid.major.x = element_line(color = "#e6e6e6", linewidth = 0.4),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(25, 25, 20, 25)
        )
      
      # Save horizontal version
      ggsave(
        filename = file.path(output_dir, "isotope_violin_plot_k6_horizontal.png"),
        plot = p_violin_horizontal,
        width = 10,
        height = 8,
        dpi = 600,
        bg = "white"
      )
      
      cat(sprintf("✓ Saved horizontal violin plot to: %s\n", 
                  file.path(output_dir, "isotope_violin_plot_k6_horizontal.png")))
    }
  }
  
  # ============================================================================
  # EXPORT RESULTS
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("EXPORTING RESULTS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Export elbow method metrics
  elbow_metrics <- data.frame(
    k = K_range,
    wss = wss,
    between_ss = between_ss,
    silhouette = silhouette_avg,
    wss_decrease = wss_change,
    wss_pct_decrease = wss_pct_change
  )
  
  write.csv(elbow_metrics, 
            file.path(output_dir, "elbow_method_metrics_combined_upper_kusko.csv"),
            row.names = FALSE)
  
  cat(sprintf("✓ Saved elbow metrics to: %s\n", 
              file.path(output_dir, "elbow_method_metrics_combined_upper_kusko.csv")))
  
  # Export cluster assignments for each k
  for (k in optimal_k_values) {
    write.csv(cluster_results[[as.character(k)]], 
              file.path(output_dir, sprintf("cluster_assignments_k%d_combined_upper_kusko.csv", k)),
              row.names = FALSE)
    
    cat(sprintf("✓ Saved k=%d cluster assignments to: %s\n", 
                k, file.path(output_dir, sprintf("cluster_assignments_k%d_combined_upper_kusko.csv", k))))
  }
  
  # Export management unit summary
  write.csv(mgmt_summary, 
            file.path(output_dir, "management_unit_isotope_summary_combined_upper_kusko.csv"),
            row.names = FALSE)
  
  cat(sprintf("✓ Saved isotope summary to: %s\n", 
              file.path(output_dir, "management_unit_isotope_summary_combined_upper_kusko.csv")))
  
  # ============================================================================
  # RECOMMENDATIONS
  # ============================================================================
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("CLUSTERING RECOMMENDATIONS\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Find the "elbow"
  wss_2nd_deriv <- diff(wss_change[-1])
  recommended_k <- which.max(wss_2nd_deriv) + 2
  
  cat(sprintf("Recommended number of clusters (based on elbow): k = %d\n", recommended_k))
  
  # Find k with highest silhouette score
  best_sil_k <- K_range[which.max(silhouette_avg)]
  cat(sprintf("Best k by silhouette score: k = %d (score = %.3f)\n", 
              best_sil_k, max(silhouette_avg, na.rm = TRUE)))
  
  # Percentage of variance explained
  cat("\nVariance Explained by Number of Clusters:\n")
  for (i in seq_along(K_range)) {
    cat(sprintf("  k = %2d: %.1f%%\n", K_range[i], variance_explained[i]))
  }
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("ANALYSIS COMPLETE\n")
  cat(rep("=", 80), "\n", sep = "")
  cat(sprintf("\nAll results saved to: %s/\n", output_dir))
  cat("\nGenerated files:\n")
  cat("  - elbow_method_analysis_combined_upper_kusko.png (4-panel visualization)\n")
  cat("  - elbow_method_metrics_combined_upper_kusko.csv (numerical results)\n")
  cat("  - clusters_k3_combined_upper_kusko.png, clusters_k4_combined_upper_kusko.png, etc.\n")
  cat("  - cluster_assignments_k3_combined_upper_kusko.csv, k4_combined_upper_kusko.csv, etc.\n")
  cat("  - basin_map_k3_combined_upper_kusko.png, basin_map_k4_combined_upper_kusko.png, etc.\n")
  cat("  - isotope_violin_plot_k6_publication.png (publication-quality violin plot)\n")
  cat("  - isotope_violin_plot_k6_annotated.png (with sample sizes)\n")
  cat("  - isotope_violin_plot_k6_horizontal.png (alternative horizontal layout)\n")
  cat("  - management_unit_isotope_summary_combined_upper_kusko.csv\n")
  cat("\nNOTE: Upper Kusko Main now includes all S. Fork Kusko data\n")
  cat("NOTE: Violin plots are 600 DPI publication-quality figures for k=6\n")
  cat("\n")
  
  # Return results invisibly
  invisible(list(
    elbow_metrics = elbow_metrics,
    cluster_results = cluster_results,
    mgmt_summary = mgmt_summary,
    recommended_k = recommended_k,
    best_sil_k = best_sil_k
  ))
  
}

################################################################################
# EXECUTE ANALYSIS
################################################################################

cat("\n")
cat("################################################################################\n")
cat("# ISOTOPE CLUSTERING ANALYSIS - USING ALL VALUES\n")
cat("# South Fork Kuskokwim COMBINED into Upper Kuskokwim\n")
cat(paste0("# Total isotope values analyzed: ALL individual predictions\n"))
cat("# Management units: 18 units (S. Fork merged into Upper Kusko Main)\n")
cat("# Exclusions: Johnson & Lower Kusko\n")
cat("################################################################################\n")

# Run the analysis
results <- run_isotope_clustering()

cat("\n")
cat("################################################################################\n")
cat("# ANALYSIS FINISHED\n")
cat("################################################################################\n")
cat("\nTo access results programmatically, assign to variable:\n")
cat("  results <- source('isotope_clustering_combined_upper_kusko.R')$value\n")
cat("\n")