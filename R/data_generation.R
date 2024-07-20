#' Generate Right Tail Data Based on Log-Normal Distribution
#'
#' This function generates data following a right-tailed distribution, specifically
#' a log-normal distribution, and then scales and truncates it to fit within a specified range.
#' The data is rounded to 2 decimal places.
#'
#' @param n The number of data points to generate. Must be a positive integer.
#' @param mean The mean of the log-normal distribution before scaling. Must be a positive number.
#' @param sd The standard deviation of the log-normal distribution before scaling. Must be a positive number.
#' @param min_val The minimum value in the scaled data range. Must be numeric and less than `max_val`.
#' @param max_val The maximum value in the scaled data range. Must be numeric and greater than `min_val`.
#'
#' @return A numeric vector of length `n` containing the generated data, scaled to the range [`min_val`, `max_val`],
#' and truncated to ensure all values fall within this range. Each value is rounded to 2 decimal places.
#'
#' @examples
#' generated_data <- generate_right_tail_data_ln(n = 100, mean = 50, sd = 10, min_val = 20, max_val = 100)
#' hist(generated_data)
#'
#' @export
generate_right_tail_data_ln <- function(n, mean, sd, min_val, max_val) {
  # Parameter validation
  if (!is.numeric(n) || n <= 0 || floor(n) != n) stop("n must be a positive integer.")
  if (!is.numeric(mean) || mean <= 0) stop("mean must be a positive number.")
  if (!is.numeric(sd) || sd <= 0) stop("sd must be a positive number.")
  if (!is.numeric(min_val) || !is.numeric(max_val) || min_val >= max_val) stop("min_val must be less than max_val and both must be numeric.")
  
  # Calculate log-normal parameters
  meanlog <- log(mean^2 / sqrt(sd^2 + mean^2))
  sdlog <- sqrt(log(1 + sd^2 / mean^2))
  
  # Generate log-normal distribution
  data <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  
  # Scale data to the specified range
  data_scaled <- (data - min(data)) / (max(data) - min(data)) * (max_val - min_val) + min_val
  
  # Round data to 2 decimal places
  data_rounded <- round(data_scaled, digits = 2)
  
  # Truncate data to ensure it falls within the specified range
  data_truncated <- pmax(pmin(data_rounded, max_val), min_val)
  
  return(data_truncated)
}

#' Generate Left Tail Data Based on Inverted Log-Normal Distribution
#'
#' This function generates data following a left-tailed distribution by inverting a log-normal distribution,
#' then scales and truncates it to fit within a specified range. The data is rounded to 2 decimal places.
#'
#' @param n The number of data points to generate. Must be a positive integer.
#' @param mean The mean of the original log-normal distribution before inversion and scaling. 
#' Must be a positive number.
#' @param sd The standard deviation of the original log-normal distribution before inversion and scaling. 
#' Must be a positive number.
#' @param min_val The minimum value in the scaled data range. Must be numeric and less than `max_val`.
#' @param max_val The maximum value in the scaled data range. Must be numeric and greater than `min_val`.
#'
#' @return A numeric vector of length `n` containing the generated data, scaled to the range [`min_val`, `max_val`],
#' and truncated to ensure all values fall within this range. Each value is rounded to 2 decimal places.
#'
#' @examples
#' generated_data <- generate_left_tail_data_ln(n = 100, mean = 50, sd = 10, min_val = 20, max_val = 100)
#' hist(generated_data)
#'
#' @export
generate_left_tail_data_ln <- function(n, mean, sd, min_val, max_val) {
  # Parameter validation
  if(n <= 0) stop("n must be greater than 0")
  if(sd <= 0) stop("sd must be greater than 0")
  if(min_val >= max_val) stop("min_val must be less than max_val")
  
  # Log-normal parameters
  meanlog <- log(mean^2 / sqrt(sd^2 + mean^2))
  sdlog <- sqrt(log(1 + sd^2 / mean^2))
  
  # Generate log-normal distribution
  data <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  
  # Invert the distribution to create a left tail
  data <- max(data) - data
  
  # Normalize to 0-1 range
  data_normalized <- (data - min(data)) / (max(data) - min(data))
  
  # Scale to specified range and round
  data_scaled <- round(data_normalized * (max_val - min_val) + min_val, digits = 2)
  
  # Ensure data is within specified range
  data_final <- pmin(pmax(data_scaled, min_val), max_val)
  
  return(data_final)
}

#' Generate Normally Distributed Data Within a Specified Interval
#'
#' This function generates a specified number of data points from a normal distribution
#' that fall within a given interval. It ensures that all generated data points are
#' within the specified minimum and maximum values.
#'
#' @param n The number of data points to generate. Must be a positive integer.
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution. Must be a positive number.
#' @param min_val The minimum value of the interval. Must be less than `max_val`.
#' @param max_val The maximum value of the interval. Must be greater than `min_val`.
#'
#' @return A numeric vector of length `n` containing the generated data points that fall within the specified interval.
#'
#' @examples
#' generated_data <- generate_normal_n_in_interval(n = 100, mean = 0, sd = 1, min_val = -2, max_val = 2)
#' hist(generated_data)
#'
#' @export
generate_normal_n_in_interval <- function(n, mean, sd, min_val, max_val) {
  # Parameter validation
  if(n <= 0) stop("n must be greater than 0")
  if(sd <= 0) stop("sd must be greater than 0")
  if(min_val >= max_val) stop("min_val must be less than max_val")
  
  # Initialize a vector to store the data points
  filtered_data <- numeric(0)
  
  # Generate data points in batches until we reach the desired number
  while (length(filtered_data) < n) {
    # Estimate the number of points to generate, considering some will be filtered out
    batch_size <- ceiling((n - length(filtered_data)) / (pnorm(max_val, mean, sd) - pnorm(min_val, mean, sd)))
    
    # Generate a batch of data points
    data_batch <- rnorm(batch_size, mean, sd)
    
    # Filter the data points to keep only those within the interval
    data_batch <- data_batch[data_batch >= min_val & data_batch <= max_val]
    
    # Add the filtered data points to the collected data
    filtered_data <- c(filtered_data, data_batch)
  }
  
  # If more points were generated than needed, truncate the vector
  filtered_data <- head(filtered_data, n)
  
  return(filtered_data)
}

#' Calculate Proportion of Data Within Standard Deviation Multiples from the Mean
#'
#' This function calculates the proportion of data points within a specified number of standard deviations from the mean.
#' It supports both single and multiple standard deviation values. For multiple values, it returns a detailed list
#' of proportions and their corresponding bounds.
#'
#' @param data A numeric vector containing the data points.
#' @param sd_multiple A numeric value or vector indicating the multiples of the standard deviation to consider
#' from the mean. Must be non-negative.
#'
#' @return If `sd_multiple` is a single number, returns a single numeric value representing the proportion of data
#' points within that range. If `sd_multiple` is a vector, returns a list where each element corresponds to a
#' `sd_multiple` value and contains the mean, standard deviation, lower bound, upper bound, and proportion of data
#' points within the range.
#'
#' @examples
#' data <- rnorm(1000)
#' single_proportion <- proportion_around_the_mean(data, 1)
#' multiple_proportions <- proportion_around_the_mean(data, c(1, 2, 3))
#'
#' @export
proportion_around_the_mean <- function(data, sd_multiple = 1) {
  # Parameter validation
  if (!is.numeric(data)) stop("Data must be a numeric vector.")
  if (any(sd_multiple < 0)) stop("sd_multiple must be non-negative.")
  
  # Calculate mean and standard deviation
  mean_val <- mean(data, na.rm = TRUE)
  sd_val <- sd(data, na.rm = TRUE)
  
  # Initialize a list to store results if sd_multiple is a vector
  results <- list()
  
  # Calculate the proportion for each sd_multiple value
  for (sd_mult in sd_multiple) {
    lower_bound <- mean_val - sd_mult * sd_val
    upper_bound <- mean_val + sd_mult * sd_val
    
    # Calculate the proportion within the range
    proportion <- mean(data >= lower_bound & data <= upper_bound, na.rm = TRUE)
    
    # Store results in a list
    results[[paste0(sd_mult, " SD")]] <- list(
      Mean = mean_val,
      StandardDeviation = sd_val,
      LowerBound = lower_bound,
      UpperBound = upper_bound,
      Proportion = proportion
    )
  }
  
  # If sd_multiple is a single number, return just the proportion; otherwise, return the list
  if (length(sd_multiple) == 1) {
    return(results[[1]]$Proportion)
  } else {
    return(results)
  }
}
