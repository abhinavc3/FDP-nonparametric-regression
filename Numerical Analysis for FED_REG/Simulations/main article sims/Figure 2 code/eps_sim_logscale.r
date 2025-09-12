# Load necessary libraries
library(wavethresh)

# Set working directory
setwd("/home/stat/lassev/PrivateNonparametricRegression")
#setwd("/Users/lassev/Dropbox (Penn)/Distributed comm. and privacy constraints/Federated Learning for Nonparametric Regression/Simulations/")

# Include wavelet helper functions
source("wavelet_helper_functions.r")

# Set seed for reproducibility
set.seed(2024)

# Simulation parameters
epsilon_values <- 10^seq(log10(0.1), log10(1000), length.out = 20)  # Epsilon ranging from 0.1 to 1000 on a log scale
machine_counts <- c(1, 100, 250, 2000)
N <- 2000  # Total number of observations
num_simulations <- 1000  # Reduced number of simulations per configuration for computational feasibility
t_point <- 0.5  # Point at which to calculate the squared error

# Wavelet and signal parameters
S <- 4  # Smoothness level of the wavelets
grid_size <- 2^10
max_level <- log2(grid_size)
wavelet_family <- "DaubExPhase"
boundary <- "interval"
s <- 2  # Smoothness of the ground truth signal

# Compute c_psi
c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size, 
  family = wavelet_family, bc = boundary, filter_number = S))

# x-values for plotting
x_values <- seq(0, 1, length.out = grid_size)

# Function to calculate index of point t
get_index_t <- function(t_point, x_values) {
  which.min(abs(x_values - t_point))
}

index_t <- get_index_t(t_point, x_values)

# Initialize a data frame to store results
results_df <- data.frame()

# Loop over each epsilon value
for (eps in epsilon_values) {
  # Print progress
  cat("Epsilon:", eps, "\n")

  # Loop over each machine count
  for (m in machine_counts) {
    # Print progress
    cat("  Machine count:", m, "\n")

    ns <- rep(N / m, m)  # Local observations per machine

    # Ensure ns are integers
    ns <- floor(ns)
    N_actual <- sum(ns)

    # Storage for MSEs
    mse_total_list <- numeric(num_simulations)
    mse_point_list <- numeric(num_simulations)

    for (sim in 1:num_simulations) {
      # Print progress every 10 simulations
      if (sim %% 10 == 0) {
        cat("    Simulation:", sim, "\n")
      }

      # Generate ground truth signal
      signal <- 10 * sqrt(grid_size) * besov_smooth_function(grid_size = grid_size, s = s, p = 2, q = 2)

      # Generate data
      X <- sample(1:grid_size, N_actual, replace = TRUE)
      Y <- signal[X] + rnorm(N_actual, mean = 0, sd = 1)

      # Privacy budgets per machine
      pb <- rep(eps, m)

      # Compute parameters
      L_max <- compute_L_max(ns, pb, s)
      tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max)

      # Estimate signal
      estimated_signal <- federated_estimator(ns, pb, s, Y, X, grid_size, max_level, 
        wavelet_family, boundary, S, tau, c_psi)

      # Calculate MSE over the entire signal
      mse_total <- mean((signal - estimated_signal)^2)
      mse_total_list[sim] <- mse_total

      # Calculate squared error at point t_point
      mse_point <- (signal[index_t] - estimated_signal[index_t])^2
      mse_point_list[sim] <- mse_point
    }

    # Calculate mean and standard deviation for total MSE
    mse_mean_total <- mean(mse_total_list)
    mse_sdev_total <- sd(mse_total_list)

    # Calculate mean and standard deviation for point MSE
    mse_mean_point <- mean(mse_point_list)
    mse_sdev_point <- sd(mse_point_list)

    # Store results in the data frame
    results_df <- rbind(results_df, data.frame(
      epsilon = eps,
      machine_count = m,
      mse_total = mse_mean_total,
      mse_total_sdev = mse_sdev_total,
      mse_point = mse_mean_point,
      mse_point_sdev = mse_sdev_point
    ))
  }
}

# Convert machine_count to a factor for plotting
#results_df$machine_count <- factor(results_df$machine_count)
results_df$machine_count <- factor(results_df$machine_count,
                                   levels = machine_counts,
                                   labels = paste0("m = ", machine_counts))

save(results_df, file = paste0(N,"N_",num_simulations,"simulation_eps_logscale.RData"))

make_plots <- FALSE
if(make_plots){
library(ggplot2)
#quartz()
load(paste0(num_simulations,"simulation_eps_logscale.RData"))

# Plotting Total Risk vs. Epsilon
# Define a vector of line types
line_types <- c("solid", "dashed", "dotted", "dotdash")

# Define a vector of custom colors
colors <- c("red", "blue", "green", "purple")

# Plotting Total Risk vs. Epsilon with enhanced distinguishability and smooth lines
ggplot(data = results_df, aes(x = epsilon, y = mse_total,
                              color = machine_count, linetype = machine_count)) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, size = 2) +  # Adds a smoothed line
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  labs(title = "",
       x = "epsilon",
       y = expression(Log[10] ~ "IMSE"),
       color = "Machine Count",
       linetype = "Machine Count") +
  theme_minimal() +
  theme(legend.position = c(1, 1),  # Move the legend to the top-right
      legend.justification = c("right", "top"),  # Align it to the top-right
    legend.title = element_blank(),
      legend.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20))

ggsave("IMSE_eps_asymptotics.pdf", width = 8, height = 6)

# Plotting Pointwise Risk vs. Epsilon with enhanced distinguishability and smooth lines
ggplot(data = results_df, aes(x = epsilon, y = mse_point,
                              color = machine_count, linetype = machine_count)) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, size = 2) +  # Adds a smoothed line
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  labs(title = "",
       x = "epsilon",
       y = expression(Log[10] ~ "Pointwise MSE at t = 0.5"),
       color = "Machine Count",
       linetype = "Machine Count") +
  theme_minimal() +
  theme(legend.position = c(1, 1),  # Move the legend to the top-right
      legend.justification = c("right", "top"),  # Align it to the top-right
    legend.title = element_blank(),
      legend.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20))

ggsave("pointwise_eps_asymptotics.pdf", width = 8, height = 6)


}
