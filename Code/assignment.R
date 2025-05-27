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
      
      # # Make a data frame with assign and the name column from edges 
      # assign_df <- data.frame(assign = assign, name = edges$Name)
      # 
      # # Find the top 5% of values 
      # top_5_percent <- quantile(assign_df$assign, 0.95, na.rm = TRUE)
      # 
      # # filter the values above the top 5% 
      # assign_df <- assign_df %>% filter(assign >= top_5_percent)
      # 
      # # find the name with the most segments in the top 5% 
      # top_name <- assign_df %>% group_by(name) %>% 
      #   summarise(count = n()) %>% 
      #   arrange(desc(count)) %>% 
      #   slice(1) %>% 
      #   pull(name)
      # 
      # hucs<- edges$Name
      # 
      # # find the indices of the top name in the edges data frame
      # top_indices <- which(hucs == top_name)
      # non_top_indices <- which(hucs != top_name)
      # 
      # assign[non_top_indices] <- 0 
      # 
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

#' Divide data into quartiles based on fixed DOY dates (UPDATED)
#'
#' @param natal_data Data frame with natal origins data
#' @return List containing quartile subsets, breaks, and labels
divide_doy_quartiles <- function(natal_data) {
  # Fixed DOY breaks based on calendar dates
  # June 11 = DOY 162, then 10-day intervals
  # DOY 162 = June 11
  # DOY 172 = June 21 (162 + 10)
  # DOY 182 = July 1 (162 + 20) 
  # DOY 192 = July 11 (162 + 30)
  
  doy_breaks <- c(162, 172, 182, 192)
  
  # Create quartile subsets based on fixed dates
  quartile_subsets <- list()
  
  # Q1: Before June 11 (DOY < 162)
  quartile_subsets[[1]] <- natal_data %>% 
    filter(DOY < doy_breaks[1])
  
  # Q2: June 11-20 (DOY 162-171)  
  quartile_subsets[[2]] <- natal_data %>% 
    filter(DOY >= doy_breaks[1] & DOY < doy_breaks[2])
  
  # Q3: June 21-30 (DOY 172-181)
  quartile_subsets[[3]] <- natal_data %>% 
    filter(DOY >= doy_breaks[2] & DOY < doy_breaks[3])
  
  # Q4: July 1+ (DOY 182+)
  quartile_subsets[[4]] <- natal_data %>% 
    filter(DOY >= doy_breaks[3])
  
  # Create descriptive labels with actual date ranges found in data
  subset_labels <- character(4)
  
  for (i in 1:4) {
    if (nrow(quartile_subsets[[i]]) > 0) {
      actual_min <- min(quartile_subsets[[i]]$DOY, na.rm = TRUE)
      actual_max <- max(quartile_subsets[[i]]$DOY, na.rm = TRUE)
      
      # Convert DOY to readable dates for labels
      min_date <- format(as.Date(actual_min, origin = "2020-01-01"), "%b %d")
      max_date <- format(as.Date(actual_max, origin = "2020-01-01"), "%b %d")
      
      subset_labels[i] <- sprintf("DOY Q%d: %d-%d (%s-%s)", 
                                  i, actual_min, actual_max, min_date, max_date)
    } else {
      # Handle empty quartiles
      expected_ranges <- c("Before Jun 11", "Jun 11-20", "Jun 21-30", "Jul 1+")
      subset_labels[i] <- sprintf("DOY Q%d: %s (No data)", i, expected_ranges[i])
    }
  }
  
  return(list(
    subsets = quartile_subsets,
    breaks = doy_breaks,
    labels = subset_labels
  ))
}