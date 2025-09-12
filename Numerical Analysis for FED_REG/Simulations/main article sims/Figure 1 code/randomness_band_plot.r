# Load necessary library
library(wavethresh)

setwd("/Users/lassev/Dropbox (Penn)/Distributed comm. and privacy constraints/Federated Learning for Nonparametric Regression/Simulations/")
source("wavelet_helper_functions.r")

# Set seed for reproducibility
set.seed(2024)

# Wavelet parameters
S <- 4  # Set the smoothness level of the wavelets
grid_size <- 2^10  # Reduced grid size for computational efficiency
max_level <- log2(grid_size)
wavelet_family <- "DaubExPhase"
boundary <- "interval"

# Calculate the maximum of the mother wavelet
c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size,
                            family = wavelet_family, bc = boundary, filter_number = S))

# Generate the true signal (Besov-smooth function)
s <- 2  # Smoothness of the ground truth
signal <- 10 * sqrt(grid_size) * besov_smooth_function(grid_size = grid_size, s = s, p = 2, q = 2, bc = boundary)

# Set up the simulation parameters
m <- 1  # Number of servers (centralized DP)
ns <- rep(500, m)  # Number of observations per server
N <- sum(ns)  # Total number of observations


# Define non-uniform probabilities over 1:grid_size
decay_rate <- 0.05
probabilities <- 0.2 + exp(-decay_rate * (1:grid_size))
probabilities <- probabilities / sum(probabilities)  # normalize to sum to 1

# Sample X with non-uniform probabilities
#X <- sample(1:grid_size, N, replace = TRUE, prob = probabilities)
# Generate data (same data for all simulations)
X <- sample(1:grid_size, N, replace = TRUE)
Y <- signal[X] + rnorm(N, mean = 0, sd = 1)

# Privacy parameters
eps_values <- c(0.5,1, 2, 5)  # Privacy budgets to consider

# Clipping parameter
tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * max_level)

# Number of simulations
num_simulations <- 100

# Prepare x-axis values
x_values <- seq(0, 1, length.out = grid_size)

# Initialize lists to store results for each eps
estimated_signals_list <- list()
mean_estimated_signals <- list()
conf_lower_list <- list()
conf_upper_list <- list()

# Colors and line types for plotting
colors <- c("red", "green", "purple", "orange")
line_types <- c(1, 2, 3, 4)

# Run simulations for each eps value
for (idx in 1:length(eps_values)) {
  eps <- eps_values[idx]
  pb <- rep(eps, m)  # Privacy budget per server
  
  # Compute L_max based on privacy and data parameters
  L_max <- compute_L_max(ns, pb, s)
  
  # Initialize a matrix to store the estimated signals for current eps
  estimated_signals <- matrix(0, nrow = grid_size, ncol = num_simulations)
  
  # Run the simulations
  for (sim in 1:num_simulations) {
    # Estimate the signal privately
    estimated_signal <- federated_estimator(ns, pb, s, Y = Y, X = X, grid_size = grid_size,
                                            max_level = max_level, wavelet_family = wavelet_family,
                                            boundary = boundary, S = S, tau = tau, c_psi = c_psi)
    # Store the estimated signal
    estimated_signals[, sim] <- estimated_signal
  }
  
  # Store the estimated signals matrix in the list
  estimated_signals_list[[idx]] <- estimated_signals
  
  # Calculate the pointwise mean of the estimated signals
  mean_estimated_signal <- rowMeans(estimated_signals)
  mean_estimated_signals[[idx]] <- mean_estimated_signal
  
  # Calculate the 95% confidence intervals (2.5th and 97.5th percentiles)
  conf_lower <- apply(estimated_signals, 1, quantile, probs = 0.025)
  conf_upper <- apply(estimated_signals, 1, quantile, probs = 0.975)
  conf_lower_list[[idx]] <- conf_lower
  conf_upper_list[[idx]] <- conf_upper
  print(idx)
}

# Compute the overall y-limits across all confidence bands and the true signal
all_conf_lower <- unlist(conf_lower_list)
all_conf_upper <- unlist(conf_upper_list)

y_min <- min(c(signal, all_conf_lower))
y_max <- max(c(signal, all_conf_upper))

# Optionally, add some padding to the y-axis limits
y_padding <- 0.1 * (y_max - y_min)
ylim_range <- c(y_min - y_padding, y_max + y_padding)

# For PDF format
pdf("randomness_over_privacy_mechanism.pdf", width = 8, height = 6)

# Plotting the results with adjusted y-limits
plot(x = x_values, y = signal, type = "l", col = "blue", lwd = 2,
     main = "",
     xlab = "X-values", ylab = "Y-values", ylim = ylim_range)
lines(x = x_values, y = signal, col = "blue", lwd = 2)

# For each eps, plot the mean estimated signal and confidence bands
for (idx in 1:length(eps_values)) {
  eps <- eps_values[idx]
  mean_estimated_signal <- mean_estimated_signals[[idx]]
  conf_lower <- conf_lower_list[[idx]]
  conf_upper <- conf_upper_list[[idx]]
  
  # Plot the mean estimated signal
  lines(x = x_values, y = mean_estimated_signal, col = colors[idx], lwd = 2, lty = line_types[idx])
  
  # Plot the confidence bands
  lines(x = x_values, y = conf_lower, col = colors[idx], lwd = 1, lty = line_types[idx])
  lines(x = x_values, y = conf_upper, col = colors[idx], lwd = 1, lty = line_types[idx])
  
  # Optional: Shade the confidence band area
  polygon(c(x_values, rev(x_values)), c(conf_upper, rev(conf_lower)),
          col = adjustcolor(colors[idx], alpha.f = 0.1), border = NA)
}

# Add legend
legend_labels <- paste("eps =", eps_values)
legend("topright", legend = c("True Signal", legend_labels),
       col = c("blue", colors), lwd = 2, lty = c(1, line_types), bg = 'white')
dev.off()

# Plot only the first curve and the scatter points
pdf("scatter_and_first_curves.pdf", width = 8, height = 6)

# Plotting the results with adjusted y-limits
plot(x = x_values, y = signal, type = "l", col = "blue", lwd = 4,
     main = "",
     xlab = "X-values", ylab = "Y-values", ylim = ylim_range)
lines(x = x_values, y = signal, col = "blue", lwd = 4)
points(x = X/grid_size, y = Y, col = "black", bg = NA, pch = 1)

# For each eps, plot the mean estimated signal and confidence bands
for (idx in 1:length(eps_values)) {
  eps <- eps_values[idx]
  first_estimated_signal <- estimated_signals_list[[idx]][,1]
  
  # Plot the mean estimated signal
  lines(x = x_values, y = first_estimated_signal, col = colors[idx], lwd = 4, lty = line_types[idx])
  
}

# Add legend
legend_labels <- paste("eps =", eps_values)
legend("topright", legend = c("True Signal", legend_labels),
       col = c("blue", colors), lwd = 2, lty = c(1, line_types), bg = 'white')
dev.off()



# Save as a high-quality PDF
pdf("Wavelet_Estimation_Privacy_Mechanism.pdf", width = 8, height = 6)

# Plot the true signal
plot(x = x_values, y = signal, type = "l", col = "black", lwd = 2, lty = 1,
     xlab = "X-values", ylab = "Y-values", ylim = ylim_range)

# Plot the mean estimated signals and confidence intervals for each eps value
for (idx in 1:length(eps_values)) {
  mean_estimated_signal <- mean_estimated_signals[[idx]]
  conf_lower <- conf_lower_list[[idx]]
  conf_upper <- conf_upper_list[[idx]]

  # Plot the mean estimated signal
  lines(x = x_values, y = mean_estimated_signal, col = "black", lwd = 2, lty = line_types[idx])

  # Plot the confidence bands
  lines(x = x_values, y = conf_lower, col = "black", lwd = 1, lty = line_types[idx])
  lines(x = x_values, y = conf_upper, col = "black", lwd = 1, lty = line_types[idx])

  # Optional: Add shaded confidence interval bands (set alpha transparency)
  polygon(c(x_values, rev(x_values)), c(conf_upper, rev(conf_lower)),
          col = adjustcolor("gray", alpha.f = 0.2), border = NA)
}

# Add a legend
legend_labels <- paste("eps =", eps_values)
legend("topright", legend = c("True Signal", legend_labels),
       col = "black", lwd = 2, lty = c(1, line_types), bty = "n")

# Close the PDF device
dev.off()

# Plot only the first curve and the scatter points in JASA style
pdf("scatter_and_first_curves.pdf", width = 8, height = 6)

# Plot the true signal with thicker line and different style (solid)
plot(x = x_values, y = signal, type = "l", col = "black", lwd = 3, lty = 1,  # True signal is thicker and solid
     xlab = "X-values", ylab = "Y-values", ylim = ylim_range)

# Add scatter points (scaled X), smaller and lighter
points(x = X/grid_size, y = Y, col = adjustcolor("black", alpha.f = 0.3), 
       pch = 16, cex = 0.5)  # Smaller points with lighter color

# For each eps, plot the first estimated signal
for (idx in 1:length(eps_values)) {
  first_estimated_signal <- estimated_signals_list[[idx]][,1]
  
  # Plot the estimated signal with varying line types (black for JASA style)
  lines(x = x_values, y = first_estimated_signal, col = "black", lwd = 2, lty = line_types[idx])
}

# Add legend with clear distinction between true signal and eps lines
legend_labels <- paste("eps =", eps_values)
legend("topright", legend = c("True Signal", legend_labels),
       col = "black", lwd = c(3, 2, 2, 2, 2), lty = c(1, line_types), bty = "n")  # True signal thicker

dev.off()  # Close the PDF device


### GGPLOT2 version ####

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Prepare the data for ggplot2
signal_df <- data.frame(x = x_values, y = signal, eps = "Non-private")
data_points_df <- data.frame(x = X / grid_size, y = Y)
first_estimates_df <- data.frame(
  x = rep(x_values, times = length(eps_values)),
  y = unlist(lapply(estimated_signals_list, function(estimates) estimates[, 1])),
  eps = as.character(rep(eps_values, each = length(x_values)))
)

# Combine signal_df and first_estimates_df
combined_df <- bind_rows(signal_df, first_estimates_df)

# Adjust the factor levels to include "Non-private" first
levels_eps <- c("Non-private", as.character(eps_values))
combined_df$eps <- factor(combined_df$eps, levels = levels_eps)

# Subsample points for better readability (e.g., every 10th point)
first_estimates_df_subsampled <- combined_df %>%
  filter(eps != "Non-private") %>%
  filter(x %in% x_values[seq(1, length(x_values), by = 10)])

data_points_df_subsampled <- data_points_df %>%
  slice(seq(1, n(), by = 10))

# Define color palette, line types, and shape types
color_palette <- c("blue", "red", "green", "purple", "orange")
line_types <- c("solid", "dashed", "dotted", "dotdash", "twodash")  # Different line types
shape_types <- c(15, 16, 17, 18, 19)  # NA shape for Non-private

# Create the plot with ggplot2
ggplot() +
  # Plot the true signal and first estimated signals
  geom_line(data = combined_df, aes(x = x, y = y, color = eps, linetype = eps), size = 1.2) +
  # Add subsampled scatter points for data points (suppress legend)
  geom_point(data = data_points_df, aes(x = x, y = y), color = "#504e4e", size = 1.5, show.legend = FALSE) +
  # Add subsampled points to each estimated signal line (excluding Non-private)
  geom_point(data = first_estimates_df_subsampled, 
             aes(x = x, y = y, color = eps, shape = eps), size = 2) +
  # Set color, line type, and shape scales
  scale_color_manual(values = color_palette, name = "Method") +
  scale_linetype_manual(values = line_types, name = "Method") +
  scale_shape_manual(values = shape_types, name = "Method", na.translate = FALSE) +
  # Adjust the guides to combine legends
  guides(color = guide_legend(override.aes = list(shape = c(NA, shape_types[-1]))),
         linetype = guide_legend(),
         shape = guide_none()) +
  # Labels and theme
  labs(x = "X-values", y = "Y-values") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Save the plot as a PDF
ggsave("first_curves.pdf", width = 8, height = 6)


# Prepare the data in a suitable format for ggplot2
signal_df <- data.frame(x = x_values, mean = signal, eps = "Non-private")
conf_df <- data.frame(
  x = rep(x_values, times = length(eps_values)),
  mean = unlist(mean_estimated_signals),
  lower = unlist(conf_lower_list),
  upper = unlist(conf_upper_list),
  eps = as.character(rep(eps_values, each = length(x_values)))
)

# Combine signal_df and conf_df
combined_df <- bind_rows(signal_df, conf_df)

# Adjust the factor levels for eps
levels_eps <- c("Non-private", as.character(eps_values))
combined_df$eps <- factor(combined_df$eps, levels = levels_eps)

# Subsample points for geom_point (excluding Non-private)
conf_df_subsampled <- combined_df %>%
  filter(eps != "Non-private") %>%
  filter(x %in% x_values[seq(1, length(x_values), by = 10)])

# Define color palette, line types, and shape types
color_palette <- c("blue", "red", "green", "purple", "orange")
line_types <- c("solid", "dashed", "dotted", "dotdash", "twodash")
shape_types <- c(NA, 16, 17, 18, 19)  # NA shape for Non-private

# Create the plot with ggplot2
ggplot() +
  # Plot the mean lines (including Non-private)
  geom_line(data = combined_df, aes(x = x, y = mean, color = eps, linetype = eps), size = 1.2) +
  # Add subsampled points on the estimated means (excluding Non-private)
  geom_point(data = conf_df_subsampled, aes(x = x, y = mean, color = eps, shape = eps), size = 2) +
  # Add confidence intervals as ribbons (excluding Non-private)
  geom_ribbon(
    data = combined_df %>% filter(eps != "Non-private"),
    aes(x = x, ymin = lower, ymax = upper, fill = eps),
    alpha = 0.1,
    show.legend = FALSE  # Exclude fill from the legend
  ) +
  # Set color, fill, shape, and linetype scales
  scale_color_manual(values = color_palette, name = "Method") +
  scale_fill_manual(values = color_palette[-1], guide = FALSE) +  # Exclude Non-private from fill and hide legend
  scale_shape_manual(values = shape_types, name = "Method", na.translate = FALSE) +
  scale_linetype_manual(values = line_types, name = "Method") +
  # Adjust guides to combine legends
  guides(
    color = guide_legend(override.aes = list(
      shape = shape_types,
      linetype = line_types
    )),
    linetype = "none",
    shape = "none"
  ) +
  # Labels and theme
  labs(x = "X-values", y = "Y-values") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),   # Increase x-axis tick label size
    axis.text.y = element_text(size = 14)
  )

# Save the plot as a PDF
ggsave("Bands.pdf", width = 8, height = 6)

