library(wavethresh)
library(ggplot2)
library(reshape2)

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

# Simulation grid
c_alpha_R_grid <- 2^seq(-6, -2, 0.25)
eps_values <- c(1, 2, 5)
n_reps <- 200

# Initialize MSE matrix
mse_matrix <- matrix(0, nrow = length(eps_values), ncol = length(c_alpha_R_grid))
colnames(mse_matrix) <- paste0("c_alpha_R=", signif(c_alpha_R_grid, 2))
rownames(mse_matrix) <- paste0("eps=", eps_values)

# Main simulation loop
for (i in seq_along(eps_values)) {
  eps <- eps_values[i]
  pb <- rep(eps, m)
  L_max <- compute_L_max(ns, pb, s)
  
  for (j in seq_along(c_alpha_R_grid)) {
    c_alpha_R <- c_alpha_R_grid[j]
    mse_reps <- numeric(n_reps)
    
    for (r in 1:n_reps) {
      X <- sample(1:grid_size, N, replace = TRUE)
      Y <- signal[X] + rnorm(N)
      
      tau <- sqrt(grid_size) * c_alpha_R + sqrt((2 * s + 1) * L_max)
      
      est <- federated_estimator(ns, pb, s, Y = Y, X = X, grid_size, max_level,
                                 wavelet_family, boundary, S, tau = tau, c_psi = c_psi)
      
      mse_reps[r] <- mean((est - signal)^2)
    }
    
    mse_matrix[i, j] <- mean(mse_reps)
    print(c(i, j))
  }
}

# Reshape results for plotting
df_mse <- melt(mse_matrix)
colnames(df_mse) <- c("Epsilon", "C_alpha_R", "MSE")
df_mse$C_alpha_R <- as.numeric(sub("c_alpha_R=", "", df_mse$C_alpha_R))
df_mse$Epsilon <- as.numeric(sub("eps=", "", df_mse$Epsilon))

# Plot
ggplot(df_mse, aes(x = C_alpha_R, y = MSE, color = factor(Epsilon))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_vline(xintercept = c_psi, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = c_psi, y = max(df_mse$MSE) * 0.8,
           label = expression("true "*C[alpha*","~R]), parse = TRUE,
           angle = 90, vjust = -0.5, hjust = 0, size = 4.5, fontface = "italic") +
  scale_x_continuous(
    name = expression(C[alpha*","~R]),
    trans = "log2"
  ) +
  scale_y_log10() +
  labs(
    title = expression("Average MSE vs. " * C[alpha*","~R]),
    y = "Average MSE",
    color = expression(epsilon)
  ) +
  theme_minimal(base_size = 14)

