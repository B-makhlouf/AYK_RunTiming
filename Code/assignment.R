# assignment.R
# Functions for Bayesian assignment in watershed analysis

library(sf)
library(dplyr)
library(here)

#' Perform Bayesian assignment for a given dataset
#'
#' @param natal_data Data frame with natal origins data
#' @param edges SF object with stream edges
#' @param watershed Character: "Kusko" or "Yukon"
#' @param priors List of prior values from setup_watershed_priors()
#' @param pid_iso Vector of isotope predictions
#' @param error Vector of error values
#' @param sensitivity_threshold Threshold for assignment filtering
#' @return Matrix of assignments for each data point
perform_assignment <- function(natal_data, edges, watershed, priors, 
                               pid_iso, error, sensitivity_threshold) {
  
  # Initialize assignment matrix
  assignment_matrix <- matrix(NA, nrow = length(pid_iso), ncol = nrow(natal_data))
  
  # Process each data point
  for (i in 1:nrow(natal_data)) {
    iso_o <- as.numeric(natal_data$natal_iso[i])
    
    if (watershed == "Kusko") {
      # Kusko assignment
      assign <- (1/sqrt(2*pi*error^2)) * exp(-1*(iso_o - pid_iso)^2/(2*error^2)) * 
        priors$pid_prior * priors$StreamOrderPrior * priors$PresencePrior 
      
      # Make a data frame with assign and the name column from edges 
      assign_df <- data.frame(assign = assign, name = edges$Name)
      
      # Find the top 5% of values 
      top_5_percent <- quantile(assign_df$assign, 0.95, na.rm = TRUE)
      
      # filter the values above the top 5% 
      assign_df <- assign_df %>% filter(assign >= top_5_percent)
      
      # find the name with the most segments in the top 5% 
      top_name <- assign_df %>% group_by(name) %>% 
        summarise(count = n()) %>% 
        arrange(desc(count)) %>% 
        slice(1) %>% 
        pull(name)
      
      hucs<- edges$Name
      
      # find the indices of the top name in the edges data frame
      top_indices <- which(hucs == top_name)
      non_top_indices <- which(hucs != top_name)
      
      assign[non_top_indices] <- 0 
      
    } else if (watershed == "Yukon") {
      # Yukon assignment with genetic priors
      gen.prior <- rep(0, length = length(pid_iso))
      gen.prior[priors$LYsites] <- as.numeric(natal_data$Lower[i])
      gen.prior[priors$MYsites] <- as.numeric(natal_data$Middle[i])
      gen.prior[priors$UYsites] <- as.numeric(natal_data$Upper[i])
      
      assign <- (1/sqrt(2*pi*error^2)) * exp(-1*(iso_o - pid_iso)^2/(2*error^2)) * 
        priors$pid_prior * priors$StreamOrderPrior * gen.prior
    }
    
    # Normalize and threshold
    assign_norm <- assign / sum(assign)
    assign_rescaled <- assign_norm / max(assign_norm)
    assign_rescaled[assign_rescaled < sensitivity_threshold] <- 0
    
    # Weight by CPUE
    assignment_matrix[,i] <- assign_rescaled * as.numeric(natal_data$COratio[i])
  }
  
  return(assignment_matrix)
}

#' Process assignment matrix to get basin-scale values
#'
#' @param assignment_matrix Matrix of assignments from perform_assignment()
#' @return List containing sum, rescaled, and normalized basin assignments
process_assignments <- function(assignment_matrix) {
  # Calculate basin-scale values
  basin_assign_sum <- apply(assignment_matrix, 1, sum, na.rm = TRUE)
  basin_assign_rescale <- basin_assign_sum / sum(basin_assign_sum, na.rm = TRUE)
  basin_assign_norm <- basin_assign_rescale / max(basin_assign_rescale, na.rm = TRUE)
  
  return(list(
    sum = basin_assign_sum,
    rescale = basin_assign_rescale,
    norm = basin_assign_norm
  ))
}

#' Divide data into quartiles based on DOY
#'
#' @param natal_data Data frame with natal origins data
#' @return List containing quartile subsets, breaks, and labels
divide_doy_quartiles <- function(natal_data) {
  # Get full DOY range
  full_doy_range <- range(natal_data$DOY, na.rm = TRUE)
  
  # Create breaks at even intervals
  doy_breaks <- seq(full_doy_range[1], full_doy_range[2], length.out = 5)
  
  # Create quartile subsets
  quartile_subsets <- list()
  for (i in 1:4) {
    quartile_subsets[[i]] <- natal_data %>% 
      filter(DOY >= doy_breaks[i] & DOY < doy_breaks[i+1])
  }
  
  # Create labels
  subset_labels <- sapply(1:4, function(i) {
    sprintf("DOY Q%d: %d-%d", i, ceiling(doy_breaks[i]), floor(doy_breaks[i+1]))
  })
  
  return(list(
    subsets = quartile_subsets,
    breaks = doy_breaks,
    labels = subset_labels
  ))
}