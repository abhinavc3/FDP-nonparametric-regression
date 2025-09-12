library(ggplot2)
library(gridExtra)
source('wavelet_helper_functions.r')
source('load_data.R')
###############################
# Define Age Limits
###############################
lower_age <- 30
upper_age <- 75

###############################
# 1. Helper Functions
###############################

# Uniformly subsample data over ppX (if needed)
# ------------------------------------------------------------------
#  subsample_uniform_ppX_ppY
#     • When use_subsample = TRUE   → original uniform-bin subsample
#     • When use_subsample = FALSE  → return the data unchanged
# ------------------------------------------------------------------
subsample_uniform_ppX_ppY <- function(ppX, ppY,
                                      num_bins,
                                      draw_per_bin = NULL,
                                      use_subsample = FALSE) {
  
  # --- 0. Fast path: do nothing -----------------------------------
  if (!use_subsample) {
    return(list(ppX = ppX, ppY = ppY))
  }
  
  # --- 1. Clip covariates to the fixed age range ------------------
  ppX_clipped <- pmin(pmax(ppX, lower_age), upper_age)
  
  # --- 2. Map to [0,1] and bin ------------------------------------
  ppX_unit <- (ppX_clipped - lower_age) / (upper_age - lower_age)
  ppX_unit_bins <- cut(ppX_unit,
                       breaks = seq(0, 1, length.out = num_bins + 1),
                       include.lowest = TRUE)
  
  # --- 3. Determine how many samples per bin ----------------------
  bin_counts          <- table(ppX_unit_bins)
  max_samples_per_bin <- min(bin_counts)
  
  if (!is.null(draw_per_bin)) {
    samples_per_bin <- min(draw_per_bin, max_samples_per_bin)
    if (draw_per_bin > max_samples_per_bin) {
      warning("draw_per_bin exceeds minimum bin count. Using ",
              samples_per_bin, " samples per bin instead.")
    }
  } else {
    samples_per_bin <- max_samples_per_bin
  }

  # --- 4. Draw equal-sized subsamples from each bin ---------------
  subsample_indices <- unlist(lapply(levels(ppX_unit_bins), function(bin) {
    bin_idx <- which(ppX_unit_bins == bin)
    sample(bin_idx, samples_per_bin)
  }))

  # --- 5. Return the subsampled vectors ---------------------------
  return(list(ppX = ppX[subsample_indices],
              ppY = ppY[subsample_indices]))
}

# Preprocess data:
# - Fix the age range to [lower_age, upper_age] and map these to (0,1) to obtain the x_values.
# - Standardize y on each server using that server's mean and standard deviation.
# - For x on each server, map to the closest grid point.
# - Return the processed data along with the original ppX (for plotting) and the server-specific y scaling parameters.
preprocess_data <- function(ppX, ppY, grid_size) {
  # Fix the range to [lower_age, upper_age]
  ppX_clipped <- pmin(pmax(ppX, lower_age), upper_age)
  
  # Standardize Y using the server-specific mean and standard deviation.
  y_mean <- mean(ppY)
  y_sd <- sd(ppY)
  
  # Map ppX to [0,1] using the fixed range [lower_age, upper_age]
  X_rescaled <- (ppX_clipped - lower_age) / (upper_age - lower_age)
  # Map to grid points: find the closest grid point for each x value.
  X_mapped <- round(X_rescaled * (grid_size - 1)) + 1
  order_index <- order(X_mapped)
  X_sorted <- X_mapped[order_index]
  Y_sorted <- ((ppY - y_mean) / y_sd)[order_index]
  ppX_sorted <- ppX_clipped[order_index]  # keep the clipped age values in sorted order
  
  # Create x_values on the unit interval mapped back to [lower_age, upper_age]
  x_values_unit <- seq(0, 1, length.out = grid_size)
  x_values <- x_values_unit * (upper_age - lower_age) + lower_age
  
  return(list(X_sorted = X_sorted,
              Y_sorted = Y_sorted,
              x_values = x_values,
              ppX_sorted = ppX_sorted,
              y_mean = y_mean,
              y_sd = y_sd))
}

# (Assumed) Estimate the wavelet signals for a given privacy level.
# [Other helper functions such as mother_wavelet, sensitivity,
# private_estimate_coeffD, private_estimate_coeffC, inverse_wavelet,
# and compute_L_max are assumed to be defined elsewhere.]

# Run the complete federated analysis pipeline for a given data source.
run_federated_analysis <- function(data, num_bins, grid_size, wavelet_family, boundary, S, s) {
  # Uniformly subsample the data for this machine/source.
  subsampled <- subsample_uniform_ppX_ppY(data$age, data$thalach, num_bins = num_bins)
  ppX_sub <- subsampled$ppX
  ppY_sub <- subsampled$ppY
  
  # Preprocess the data.
  preproc <- preprocess_data(ppX_sub, ppY_sub, grid_size)
  
  return(preproc)
}

###############################
# 2. Federated Analysis on Combined Data
###############################

# Parameters
num_bins <- 1        # for uniform subsampling per source
grid_size <- 1024    # power-of-2 grid
wavelet_family <- "DaubExPhase"
boundary <- "interval"
S <- 4               # filter number
s <- 3               # smoothness parameter

# Define the privacy budgets to be evaluated.
privacy_values <- c(Inf, 0.5, 1, 2, 4)

# Split combined data by source (each source = one machine)
data_list <- split(data_filtered, data_filtered$source)

# Preprocess each source individually.
preproc_list <- lapply(data_list, function(df) {
  run_federated_analysis(df, num_bins, grid_size, wavelet_family, boundary, S, s)
})

# Get the number of observations per machine.
ns <- sapply(preproc_list, function(res) length(res$X_sorted))

# Combine preprocessed data from all sources.
X_fed <- unlist(lapply(preproc_list, function(res) res$X_sorted))
Y_fed <- unlist(lapply(preproc_list, function(res) res$Y_sorted))
# For plotting, combine the original age values from all sources.
ppX_combined <- unlist(lapply(preproc_list, function(res) res$ppX_sorted))

# Compute global x_values using the fixed age range [lower_age, upper_age]
x_values <- seq(lower_age, upper_age, length.out = grid_size)

# Additional federated parameters (these remain as in the original code).
max_level <- log2(grid_size)
c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size,
                            family = wavelet_family, bc = boundary, filter_number = S))
L_max <- compute_L_max(ns, rep(10, length(ns)), s)  # initial default pb (will be overwritten in loop)
tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max)

# Loop through each privacy budget, supply it via pb, run the federated estimator, and store the results.
# The estimated signal is then rescaled back to the original scale for each server separately.
federated_estimates <- list()
for (pb_val in privacy_values) {
  # Create a privacy vector for all machines.
  pb_vector <- rep(pb_val, length(ns))
  
  # For an infinite privacy budget, set tau to Inf; otherwise, use the computed tau.
  current_tau <- if (is.infinite(pb_val)) Inf else tau
  
  # Call the federated estimator (assumed to accept pb via the pb argument).
  est_signal <- federated_estimator(ns, pb_vector, s, Y = Y_fed, X = X_fed, grid_size = grid_size,
                                    max_level = max_level, wavelet_family = wavelet_family,
                                    boundary = boundary, S = S, tau = current_tau, c_psi = c_psi)
  
  # Split the estimated signal back to the original servers based on the number of observations.
  # Here we assume that the estimator returns a vector of length equal to grid_size.
  # For each server, rescale using the server-specific y scaling parameters.
  est_signal_list <- list()
  for (i in seq_along(preproc_list)) {
    server_y_mean <- preproc_list[[i]]$y_mean
    server_y_sd <- preproc_list[[i]]$y_sd
    # Rescale the estimated signal from the federated estimator back to the original scale.
    est_signal_list[[i]] <- est_signal * server_y_sd + server_y_mean
  }
  
  # Use a key for the privacy budget.
  key <- if (is.infinite(pb_val)) "Inf" else as.character(pb_val)
  federated_estimates[[key]] <- est_signal_list
}

###############################
# 3. Plot Federated Estimator Results using ggplot2
###############################
compute_R2 <- function(predicted, observed) {
  na_count <- sum(is.na(predicted)) + sum(is.na(observed))
  cat("Number of NA points:", na_count, "\n")
  
  rss <- sum((observed - predicted)^2, na.rm = TRUE)
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  return(1 - rss/tss)
}



# Function to plot estimated signal with true data
plot_estimated_signal_with_true_data <- function(x_values, true_x, true_y, federated_estimates, title) {
  # Create a data frame for the true data points.
  true_data <- data.frame(x = pmin(pmax(true_x, lower_age), upper_age),
                          y = true_y,
                          stringsAsFactors = FALSE)
  
  # For each privacy level, average the estimated signals across servers.
  estimate_df <- do.call(rbind, lapply(names(federated_estimates), function(key) {
    server_estimates <- federated_estimates[[key]]
    # Average over servers at each grid point.
    avg_estimate <- rowMeans(do.call(cbind, server_estimates))
    data.frame(x = x_values,
               y = avg_estimate,
               privacy = as.character(key),
               stringsAsFactors = FALSE)
  }))
  
  # Build the ggplot.
  p <- ggplot() +
    geom_point(data = true_data, aes(x = x, y = y), color = "black", size = 1, alpha = 0.5) +
    geom_line(data = estimate_df, aes(x = x, y = y, color = privacy), size = 1) +
    labs(
         #x = paste0("Age (", lower_age, "-", upper_age, ")"),
         x = "Age",
         y = "Thalach",
         color =  expression(epsilon))+
    theme_minimal()
  
  return(list(plot = p, estimate_df = estimate_df))
}

# Generate the plot and get the estimate data frame.
plot_result <- plot_estimated_signal_with_true_data(x_values, 
                                                    data_filtered$age, 
                                                    data_filtered$thalach,
                                                    federated_estimates,
                                                    title = "Estimated Signal vs. True Data")
p_true_vs_est <- plot_result$plot
estimate_df <- plot_result$estimate_df
print(p_true_vs_est)

# Compute R^2 for each privacy level.
# We will interpolate the averaged estimated signal at the locations of the true data points.
R2_results <- lapply(unique(estimate_df$privacy), function(pb_key) {
  # Subset the estimated signal for this privacy level.
  est_sub <- subset(estimate_df, privacy == pb_key)
  # Interpolate the estimated signal at the true data x-values.
  predicted_at_true <- approx(x = est_sub$x, y = est_sub$y, xout = data_filtered$age)$y
  # Compute R^2 between the interpolated predictions and the true y-values.
  R2_value <- compute_R2(predicted_at_true, data_filtered$thalach)
  return(data.frame(Privacy = pb_key, R2 = R2_value, stringsAsFactors = FALSE))
})

# Combine R2 results into one table and print.
global_R2_table <- do.call(rbind, R2_results)
print(global_R2_table)





###############################
# Additional Code: Group-wise R^2 Evaluation
###############################

library(ggplot2)
library(gridExtra)

# Define the compute_R2 function.
compute_R2 <- function(predicted, observed) {
  na_count <- sum(is.na(predicted)) + sum(is.na(observed))
  cat("Number of NA points:", na_count, "\n")
  
  rss <- sum((observed - predicted)^2, na.rm = TRUE)
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  return(1 - rss/tss)
}

# Function to create a plot for one source and compute R^2.
plot_source_estimates <- function(source_index, x_values, preproc_list, federated_estimates) {
  # Extract the preprocessed data for the source.
  source_data <- preproc_list[[source_index]]
  true_x <- source_data$ppX_sorted
  true_y <- source_data$Y_sorted * source_data$y_sd + source_data$y_mean
  
  # Create a data frame for the true data points.
  true_df <- data.frame(x = true_x, y = true_y, stringsAsFactors = FALSE)
  
  # For each privacy level, extract the estimated signal for this source.
  est_df <- do.call(rbind, lapply(names(federated_estimates), function(key) {
    est_signal <- federated_estimates[[key]][[source_index]]
    data.frame(x = x_values,
               y = est_signal,
               privacy = as.character(key),
               stringsAsFactors = FALSE)
  }))
  
  # Compute R^2 for each privacy level.
  r2_results <- lapply(unique(est_df$privacy), function(pb_key) {
    est_sub <- subset(est_df, privacy == pb_key)
    # Interpolate the estimated signal at the true x-values.
    predicted_at_true <- approx(x = est_sub$x, y = est_sub$y, xout = true_x)$y
    r2_val <- compute_R2(predicted = predicted_at_true, observed = true_y)
    return(data.frame(Privacy = pb_key, R2 = r2_val, stringsAsFactors = FALSE))
  })
  r2_df <- do.call(rbind, r2_results)
  
  # Print R^2 values for this source.
  cat("Source", source_index, "R^2 values:\n")
  print(r2_df)
  
  # Create a subtitle with R^2 information.
  subtitle_text <- paste(apply(r2_df, 1, function(row) {
    paste0(row["Privacy"], ": ", round(as.numeric(row["R2"]), 3))
  }), collapse = ", ")
  
  # Build the ggplot.
  p <- ggplot() +
    geom_point(data = true_df, aes(x = x, y = y), color = "black", size = 1, alpha = 0.6) +
    geom_line(data = est_df, aes(x = x, y = y, color = privacy), size = 1) +
    labs(title = paste0("Source ", source_index, " Estimates"),
         subtitle = subtitle_text,
         #x = paste0("Age (", lower_age, "-", upper_age, ")"),
         x = "Age",
         y = "Thalach (Original Scale)",
         color = "Privacy Budget") +
    theme_minimal()
  
  return(list(plot = p, R2 = r2_df))
}

# Create a list of plots (and R^2 data) for each source.
num_sources <- length(preproc_list)
source_plots <- lapply(seq_len(num_sources), function(i) {
  plot_source_estimates(i, x_values, preproc_list, federated_estimates)$plot
})

# Arrange the four plots in a 2x2 grid.
grid.arrange(grobs = source_plots, ncol = 2)











