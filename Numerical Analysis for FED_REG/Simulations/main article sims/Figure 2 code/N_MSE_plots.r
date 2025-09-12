# Load necessary library
library(wavethresh)

# Set working directory and source helper functions
# Update the path as per your local environment
#setwd("/Users/lassev/Dropbox (Penn)/Distributed comm. and privacy constraints/Federated Learning for Nonparametric Regression/Simulations/")
setwd("/home/stat/lassev/PrivateNonparametricRegression")
source("wavelet_helper_functions.r")

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

# Specify the point at which to compute the squared error
t_point <- 0.5  # Point in [0,1]
# Find the index corresponding to t_point
t_index <- which.min(abs(x_values - t_point))

# Total number of observations N ranging from 10 to 1000
N_values <- seq(10, 1000, by = 10)

# Define the setups
setups <- list(
  list(label = "m=1", m_func = function(N) {1}),
  list(label = "m=round(N/10)", m_func = function(N) {max(1, floor(N/10))}),
  list(label = "m=N", m_func = function(N) {N})
)

# Initialize lists to store results for each setup
results_list <- vector("list", length(setups))
names(results_list) <- sapply(setups, function(x) x$label)

# Run simulations for each setup and each N
for (setup_idx in 1:length(setups)) {
  setup <- setups[[setup_idx]]
  label <- setup$label
  m_func <- setup$m_func
  
  # Initialize vectors to store average errors for each N
  average_total_MSEs <- numeric(length(N_values))
  average_point_MSEs <- numeric(length(N_values))
  
  # Loop over N_values
  for (N_idx in 1:length(N_values)) {
    N <- N_values[N_idx]
    m <- m_func(N)
    m <- max(1, m)  # Ensure m is at least 1
    
    # Determine the number of observations per server
    ns_per_machine <- floor(N / m)
    ns <- rep(ns_per_machine, m)
    remainder <- N - sum(ns)
    if (remainder > 0) {
      # Distribute the remainder among the first few machines
      ns[1:remainder] <- ns[1:remainder] + 1
    }
    N_actual <- sum(ns)  # Adjusted total N
    
    # Compute L_max based on privacy and data parameters
    pb <- rep(eps, m)  # Privacy budget per server
    L_max <- compute_L_max(ns, pb, s)
    
    # Initialize vectors to store MSEs for current N
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
      estimated_signal <- federated_estimator(ns, pb, s, Y = Y, X = X, grid_size = grid_size,
                                              max_level = max_level, wavelet_family = wavelet_family,
                                              boundary = boundary, S = S, tau = tau, c_psi = c_psi)
      
      # Compute the total MSE for this simulation
      total_MSEs[sim] <- mean((signal - estimated_signal)^2)
      
      s_point <- s - 0.5  # Adjust the smoothness for the point error calculation

      # Estimate the signal privately
      estimated_signal <- federated_estimator(ns, pb, s_point, Y = Y, X = X, grid_size = grid_size,
                                              max_level = max_level, wavelet_family = wavelet_family,
                                              boundary = boundary, S = S, tau = tau, c_psi = c_psi)

      # Compute the squared error at the specified point
      point_MSEs[sim] <- (signal[t_index] - estimated_signal[t_index])^2
    }
    
    # Compute the average errors for current N
    average_total_MSEs[N_idx] <- mean(total_MSEs)
    average_point_MSEs[N_idx] <- mean(point_MSEs)
  }
  
  # Store the results in the list
  results_list[[label]] <- list(
    N_values = N_values,
    average_total_MSEs = average_total_MSEs,
    average_point_MSEs = average_point_MSEs
  )
}

# Save the results to a file for later retrieval
#save(results_list, file = "simulation_results_N.RData")

# Load the results (if needed)
# load("simulation_results_N.RData")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming `results_list` is already loaded or created from your simulation code

# Prepare data for Average Total MSE
total_MSE_data <- lapply(names(results_list), function(setup_label) {
  data.frame(
    N_values = results_list[[setup_label]]$N_values,
    MSE = results_list[[setup_label]]$average_total_MSEs,
    Setup = setup_label
  )
}) %>% bind_rows()

# Prepare data for Average Point MSE
point_MSE_data <- lapply(names(results_list), function(setup_label) {
  data.frame(
    N_values = results_list[[setup_label]]$N_values,
    MSE = results_list[[setup_label]]$average_point_MSEs,
    Setup = setup_label
  )
}) %>% bind_rows()


### altnerative

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare data for Average Total MSE
total_MSE_data <- lapply(names(results_list), function(setup_label) {
  data.frame(
    N_values = results_list[[setup_label]]$N_values,
    MSE = results_list[[setup_label]]$average_total_MSEs,
    Setup = setup_label
  )
}) %>% bind_rows()

# Prepare data for Average Point MSE
point_MSE_data <- lapply(names(results_list), function(setup_label) {
  data.frame(
    N_values = results_list[[setup_label]]$N_values,
    MSE = results_list[[setup_label]]$average_point_MSEs,
    Setup = setup_label
  )
}) %>% bind_rows()

# Filter data for N ≥ 100
filtered_total_MSE_data <- total_MSE_data %>% filter(N_values >= 100)
filtered_point_MSE_data <- point_MSE_data %>% filter(N_values >= 100)

# Filter data for N ≥ 100
filtered_point_MSE_data <- point_MSE_data %>% filter(N_values >= 100)

# Add log-transformed MSE to the filtered data
filtered_point_MSE_data <- filtered_point_MSE_data %>%
  mutate(Log_MSE = log10(MSE))

# Add log-transformed MSE to the filtered data
filtered_total_MSE_data <- filtered_total_MSE_data %>%
  mutate(Log_MSE = log10(MSE))

# # Add log-transformed MSE to the filtered data
# filtered_total_MSE_data <- filtered_total_MSE_data %>% mutate(Log_MSE = log10(MSE))
# filtered_point_MSE_data <- filtered_point_MSE_data %>% mutate(Log_MSE = log10(MSE))



# Plot the log-transformed data
ggplot(filtered_total_MSE_data, aes(x = N_values, y = Log_MSE, color = Setup, linetype = Setup, shape = Setup)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Total Number of Observations (N)",
    y = expression(Log[10] ~ "of IMSE"),
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("green", "blue", "red")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 17, 18))

# Save the plot as a PDF
ggsave("average_total_MSE_vs_N_log_filtered.pdf", width = 8, height = 6)

# Plot the log-transformed data
ggplot(filtered_point_MSE_data, aes(x = N_values, y = Log_MSE, color = Setup, linetype = Setup, shape = Setup)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Total Number of Observations (N)",
    y = expression(Log[10] ~ "of pointwise MSE at t=0.5"),
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  scale_color_manual(values = c("green", "blue", "red")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 17, 18))

# Save the plot as a PDF
ggsave("average_point_MSE_vs_N_log_filtered.pdf", width = 8, height = 6)


# Plot Average Total MSE
ggplot(filtered_total_MSE_data, aes(x = N_values, y = Log_MSE, color = Setup, linetype = Setup, shape = Setup)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Total Number of Observations (N)",
    y = expression(Log[10] ~ "of IMSE"),
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  scale_color_manual(values = c("blue", "red", "green")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 17, 18))

ggsave("average_total_MSE_vs_N_log_filtered.pdf", width = 8, height = 6)

# Plot Average Point MSE
ggplot(filtered_point_MSE_data, aes(x = N_values, y = Log_MSE, color = Setup, linetype = Setup, shape = Setup)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Total Number of Observations (N)",
    y = expression(Log[10] ~ "of pointwise MSE at t=0.5"),
    title = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  scale_color_manual(values = c("blue", "red", "green")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 17, 18))

ggsave("average_point_MSE_vs_N_log_filtered.pdf", width = 8, height = 6)



