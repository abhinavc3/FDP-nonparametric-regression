
setwd("/home/stat/lassev/PrivateNonparametricRegression")
source("wavelet_helper_functions.r")

# Load necessary library
library(wavethresh)

# Set seed for reproducibility
set.seed(2024)

# Wavelet parameters
S <- 4  # Smoothness level of the wavelets
grid_size <- 2^10  # Adjusted for computational efficiency
max_level <- log2(grid_size)
wavelet_family <- "DaubExPhase"
boundary <- "interval"

# Calculate the maximum of the mother wavelet
c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size,
                            family = wavelet_family, bc = boundary, filter_number = S))

# Set up the simulation parameters
s <- 2  # Smoothness of the ground truth

# Prepare x-axis values
x_values <- seq(0, 1, length.out = grid_size)

# Privacy parameter (fixed)
eps <- 1  # Privacy budget per server

# Clipping parameter
tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * max_level)

# Number of simulations
num_simulations <- 1000  # Adjust for computational efficiency

# Total number of observations N (fixed)
N <- 2^11  # Adjust as needed

# Define n values that double each time, ranging from 1 to N (logarithmic scale)
n_values <- 2^(0:floor(log2(N)))  # This gives 1, 2, 4, ..., up to the largest power of 2 ≤ N

# Initialize vectors to store average errors for each n value
average_total_MSEs <- numeric(length(n_values))
average_point_MSEs <- numeric(length(n_values))

# Specify the point at which to compute the squared error
t_point <- 0.5  # Point in [0,1]
# Find the index corresponding to t_point
t_index <- which.min(abs(x_values - t_point))

# Run simulations for each n value
for (idx in 1:length(n_values)) {
  n <- n_values[idx]
  m <- max(1, round(N / n))  # Compute m as N / n
  ns <- rep(n, m)  # Allocate n observations per machine
  remainder <- N - sum(ns)
  if (remainder > 0) {
    ns[1:remainder] <- ns[1:remainder] + 1
  }
  N_actual <- sum(ns)  # Adjusted total N
  
  # Compute L_max based on privacy and data parameters
  pb <- rep(eps, m)  # Privacy budget per server
  L_max <- compute_L_max(ns, pb, s)
  
  # Initialize vectors to store MSEs for current n value
  total_MSEs <- numeric(num_simulations)
  point_MSEs <- numeric(num_simulations)
  
  # Run the simulations
  for (sim in 1:num_simulations) {
    # Generate the true signal
    signal <- 10 * sqrt(grid_size) * besov_smooth_function(grid_size = grid_size, s = s, p = 2, q = 2, bc = boundary)
    
    # Generate data
    X <- sample(1:grid_size, N_actual, replace = TRUE)
    Y <- signal[X] + rnorm(N_actual, mean = 0, sd = 1)
    
    # Estimate the signal privately
    estimated_signal <- federated_estimator(ns, pb, s, Y = Y, X = X, grid_size = grid_size, max_level = max_level, wavelet_family = wavelet_family, boundary = boundary, S = S, tau = tau, c_psi = c_psi)

    # Compute the total MSE for this simulation
    total_MSEs[sim] <- mean((signal - estimated_signal)^2)
    
    # Compute the squared error at the specified point
    point_MSEs[sim] <- (signal[t_index] - estimated_signal[t_index])^2
  }
  
  # Compute the average errors for current n value
  average_total_MSEs[idx] <- mean(total_MSEs)
  average_point_MSEs[idx] <- mean(point_MSEs)
}

# Save the results to an RData file
save(n_values, average_total_MSEs, average_point_MSEs, file = "n_over_m_results_log.RData")

plots <- FALSE
if(plots){
# Load the results
load("n_over_m_results_log.RData")

# Convert n_values and corresponding m values to their ratios n/m for plotting
nm_ratios <- n_values / N

# Plotting Average Total MSE vs n/m on a log-log scale
plot(nm_ratios, average_total_MSEs, log = "xy", type = "b", pch = 19, col = "blue",
     xlab = "n/m (Degree of Distribution)", ylab = "IMSE",
     main = "Average IMSE vs Degree of Distribution (log-log scale)")

# Plotting Average Point MSE vs n/m on a log-log scale
plot(nm_ratios, average_point_MSEs, log = "xy", type = "b", pch = 19, col = "red",
     xlab = "n/m (Degree of Distribution)", ylab = paste("Average Squared Error at t =", t_point),
     main = paste("Average Squared Error at t =", t_point, "vs Degree of Distribution (log-log scale)"))
}

# Saving the plot
# Saving the plot
plots <- TRUE
if (plots) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Load the results
  load("n_over_m_results_log.RData")

  # Total number of observations N (fixed)
  N_total <- 2^11  # As defined in your simulation code

  # Convert n_values and corresponding m values to their ratios n/N for plotting
  nm_ratios <- n_values / N_total  # n / N

  # Prepare the data for plotting
  plot_data <- data.frame(
    nm_ratio = nm_ratios,
    Average_Total_MSE = average_total_MSEs,
    Average_Point_MSE = average_point_MSEs
  )
  plot_data <- plot_data[-1,]

  # Reshape the data to long format for ggplot2
  plot_data_long <- plot_data %>%
    pivot_longer(cols = c(Average_Total_MSE, Average_Point_MSE),
                 names_to = "Metric",
                 values_to = "MSE")

  # Create the combined plot using ggplot2
  ggplot(plot_data_long, aes(x = nm_ratio, y = MSE, color = Metric, shape = Metric, linetype = Metric)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c("blue", "red"), labels = c("IMSE", paste("Pointwise MSE at t =", t_point))) +
    scale_shape_manual(values = c(16, 17), labels = c("IMSE", paste("Pointwise MSE at t =", t_point))) +
    scale_linetype_manual(values = c("solid", "dashed"), labels = c("IMSE", paste("Pointwise MSE at t =", t_point))) +
    labs(
      x = expression(n / N ~ "Degree of Distribution (log-scale)"),
      y = expression(Log[10] ~ "of MSE"),
      title = "",
      color = "Metric",
      shape = "Metric",
      linetype = "Metric"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20)
    )

  # Save the plot as a PDF
  ggsave("Combined_MSE_vs_nm_ratio.pdf", width = 8, height = 6)
}


