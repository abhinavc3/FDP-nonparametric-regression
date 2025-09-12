library(wavethresh)


# Set working directory
setwd("/Users/abhinav/Dropbox (Personal)/Work Stuff/Distributed comm. and privacy constraints/Federated Learning for Nonparametric Regression/Simulations/")

# Include wavelet helper functions
source("wavelet_helper_functions.r")
set.seed(2024)

# Fixed parameters
m <- 1
ns <- rep(500, m)
N <- sum(ns)
grid_size <- 2^10
s <- 1
S <- 6
max_level <- log2(grid_size)
wavelet_family <- "DaubExPhase"
boundary <- "interval"

# Generate true signal
signal <- 10 * sqrt(grid_size) * besov_smooth_function(grid_size = grid_size, s = s, p = 2, q = 2)

# Fixed wavelet constant
c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size,
                            family = wavelet_family, bc = boundary, filter_number = S))

# Sim grid
sigma_grid <- 2^seq(-1, 2, 0.25)
eps_values <- c(1, 2, 5)
n_reps <- 100

# Storage
mse_matrix <- matrix(0, nrow = length(eps_values), ncol = length(sigma_grid))
colnames(mse_matrix) <- paste0("sigma=", signif(sigma_grid, 2))
rownames(mse_matrix) <- paste0("eps=", eps_values)

# Loop over eps and sigma
for (i in seq_along(eps_values)) {
  eps <- eps_values[i]
  pb <- rep(eps, m)
  L_max <- compute_L_max(ns, pb, s)
  
  for (j in seq_along(sigma_grid)) {
    sigma <- sigma_grid[j]
    mse_reps <- numeric(n_reps)
    
    for (r in 1:n_reps) {
      # Generate new data each repetition
      X <- sample(1:grid_size, N, replace = TRUE)
      Y <- signal[X] + rnorm(N)
      
      tau <- sqrt(grid_size) * c_psi + sqrt((2*s + 1) * L_max) * sigma
      
      # Estimate
      est <- federated_estimator(ns, pb, s, Y = Y, X = X, grid_size, max_level,
                                 wavelet_family, boundary, S, tau = tau, c_psi = c_psi)
      
      # Compute MSE
      mse_reps[r] <- mean((est - signal)^2)
    }
    
    # Average over repetitions
    mse_matrix[i, j] <- mean(mse_reps)
    print(c(i,j))
  }
}

library(ggplot2)
library(reshape2)

# Assuming df_mse has these columns:
# Epsilon, Sigma, MSE

df_mse <- melt(mse_matrix)
colnames(df_mse) <- c("Epsilon", "Sigma", "MSE")
df_mse$Sigma <- as.numeric(sub("sigma=", "", df_mse$Sigma))
df_mse$Epsilon <- as.numeric(sub("eps=", "", df_mse$Epsilon))

df_mse$log2_sigma <- log2(df_mse$Sigma)


ggplot(df_mse, aes(x = Sigma, y = MSE, color = factor(Epsilon))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 1, y = max(df_mse$MSE) * 0.9, label = expression("true " * sigma), 
           angle = 90, vjust = -0.5, hjust = 0, size = 4.5, fontface = "italic") +
  scale_x_continuous(
    name = expression(sigma),
    trans = "log2"
  ) +
  labs(
    title = expression("Average MSE vs. " * sigma),
    y = "Average MSE",
    color = expression(epsilon)
  ) +
  theme_minimal(base_size = 14)

## our estimator for sigma

estimate_sigma_dyadic <- function(Y, m, eps, r =2) {
  N <- length(Y)
  stopifnot(N %% m == 0)  # Equal samples per server
  
  n <- N / m
  j_vals <- -10:10
  intervals <- lapply(j_vals, function(j) c(r^j, r^(j + 1)))
  
  # Split Y into m servers
  server_Ys <- split(Y, rep(1:m, each = n))
  
  # Store noisy proportions for each server and interval
  noisy_props <- matrix(0, nrow = m, ncol = length(j_vals))
  
  for (i in 1:m) {
    Y_i <- server_Ys[[i]]
    
    for (j in seq_along(j_vals)) {
      bounds <- intervals[[j]]
      prop <- mean(Y_i > bounds[1] & Y_i <= bounds[2])
      
      # Add Laplace noise for DP
      noisy_props[i, j] <- prop + rlaplace(1, scale = 1 / (n * eps))
    }
  }
  
  # Average across servers
  avg_noisy_props <- colMeans(noisy_props)
  
  # Select the interval with max average noisy proportion
  j_star_index <- which.max(avg_noisy_props)
  j_star <- j_vals[j_star_index]
  
  # Estimate sigma
  #sigma_hat <-  r^j_star
  sigma_hat <- 5 * r^j_star
  return(sigma_hat)
}

estimate_sigma_dyadic(Y, m =1, eps = 2, r = 1.1)
