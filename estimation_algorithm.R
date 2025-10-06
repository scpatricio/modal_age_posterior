# Load the required library for multivariate normal distribution functions
library(mvtnorm)

# Function to estimate mode distribution using Gaussian approximations
mode_dist_gauss <- function(ex, dx, x) {
  # Total exposure for the interval
  # ex: Exposure vector for each age interval
  # dx: Number of deaths in each age interval
  # x: Age intervals
  
  # Define parameters
  n <- sum(dx)              # Total number of deaths across all intervals
  K <- length(dx)           # Number of age intervals
  
  # Check if the total number of deaths (n) is zero to avoid division errors
  if (n == 0) {
    stop("Total death count (n) is zero, unable to compute probabilities.")
  }
  
  # Step 1: Calculate the force of mortality for each interval
  mu_x <- dx / ex  # Instantaneous mortality rate (force of mortality)
  
  # Step 2: Compute average fraction of interval lived by those who die (a_x)
  ax <- 1 + 1 / mu_x - 1 / (1 - exp(-1 * mu_x)) 
  
  # Step 3: Compute the probability of dying within each interval (n_q_x)
  nqx <- (1 * mu_x) / (1 + (1 - ax) * mu_x)
  nqx[K] <- 1  # Ensure the probability of dying in the last interval is 1
  
  # Step 4: Compute the probability of surviving each interval (n_p_x)
  npx <- 1 - nqx
  
  # Step 5: Compute the proportion of survivors at the beginning of each interval (l_x)
  lx <- cumprod(c(1, npx[-K]))  # Survivors as a cumulative product of survival probabilities
  
  # Step 6: Compute the proportion of deaths within each interval (d_x)
  ndx <- nqx * lx
  
  # Step 7: Calculate the probabilities for each age category
  p <- ndx  # Proportion of deaths in each interval, treated as probabilities
  
  # Step 8: Define the mean vector and covariance matrix for X
  mu_X <- n * p  # Mean vector for the multivariate normal distribution
  Sigma_X <- n * (diag(p) - p %*% t(p))  # Covariance matrix
  
  # Step 9: Initialize storage for the probability of each category being the maximum
  max_probabilities <- numeric(K)
  
  # Helper function to make covariance matrices invertible
  make_invertible <- function(matrix, epsilon = 1e-16) {
    diag(matrix) <- diag(matrix) + epsilon  # Add small value to diagonal
    return(matrix)
  }
  
  # Step 10: Loop over each category to calculate the probability that X_k is the maximum
  for (k in 1:K) {
    # Create a difference matrix for the current category
    diff_matrix <- matrix(0, nrow = K - 1, ncol = K)
    idx <- 1
    for (j in 1:K) {
      if (j != k) {
        diff_matrix[idx, k] <- 1  # 1 for the category of interest
        diff_matrix[idx, j] <- -1  # -1 for all other categories
        idx <- idx + 1
      }
    }
    
    # Covariance matrix of differences
    cov_matrix <- diff_matrix %*% Sigma_X %*% t(diff_matrix)
    
    # Ensure the covariance matrix is invertible
    if (det(cov_matrix) == 0) {
      cov_matrix <- make_invertible(cov_matrix)
    }
    
    # Mean vector and bounds for the multivariate normal distribution
    mean_vec <- as.vector(diff_matrix %*% mu_X)
    lower_bound <- rep(0, K - 1)  # Lower bound for probabilities
    
    # Calculate the probability that the current category is the maximum
    max_probabilities[k] <- mvtnorm::pmvnorm(lower = lower_bound, mean = mean_vec, sigma = cov_matrix)[1]
  }
  
  # Step 11: Normalize the probabilities to ensure they form a valid distribution
  if (sum(max_probabilities) == 0) {
    warning("Sum of max probabilities is zero; setting all probabilities to zero.")
    probabilities <- rep(0, K)
  } else {
    probabilities <- max_probabilities / sum(max_probabilities)
  }
  
  # Step 12: Create a data frame with results
  df_mode <- data.frame(
    mode = x,  # Age categories
    prob = probabilities,  # Probability for each category
    cum_prob = cumsum(probabilities)  # Cumulative probabilities
  )
  
  # Step 13: Add quantiles to the data frame
  df_mode <- data.frame(
    mode = x,
    prob = probabilities
  )
  
  # Return the data frame
  return(df_mode)
}
