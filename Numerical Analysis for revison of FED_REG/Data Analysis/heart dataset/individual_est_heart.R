library(ggplot2)
library(gridExtra)
source('wavelet_helper_functions.r')
source('load_data.R')
###############################
# 1. Helper Functions
###############################

# Uniformly subsample the data
subsample_uniform_ppX_ppY <- function(ppX, ppY, num_bins, draw_per_bin = NULL) {
  # Normalize ppX to [0,1]
  ppX_unit <- (ppX - min(ppX)) / (max(ppX) - min(ppX))
  # Create bins
  ppX_unit_bins <- cut(ppX_unit, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)
  bin_counts <- table(ppX_unit_bins)
  max_samples_per_bin <- min(bin_counts)
  
  if (!is.null(draw_per_bin)) {
    samples_per_bin <- min(draw_per_bin, max_samples_per_bin)
    if (draw_per_bin > max_samples_per_bin) {
      warning("draw_per_bin exceeds the minimum bin count. Using ", samples_per_bin, " samples per bin instead.")
    }
  } else {
    samples_per_bin <- max_samples_per_bin
  }
  
  subsample_indices <- c()
  for (bin in levels(ppX_unit_bins)) {
    bin_indices <- which(ppX_unit_bins == bin)
    sampled_indices <- sample(bin_indices, samples_per_bin)
    subsample_indices <- c(subsample_indices, sampled_indices)
  }
  
  return(list(ppX = ppX[subsample_indices], ppY = ppY[subsample_indices]))
}

# Preprocess data: rescale X to [0,1], standardize Y, map X to grid, and sort.
preprocess_data <- function(ppX, ppY, grid_size) {
  # Save original Y stats for inverse transformation
  y_mean <- mean(ppY)
  y_sd <- sd(ppY)
  
  # Rescale X to [0,1]
  X_rescaled <- (ppX - min(ppX)) / (max(ppX) - min(ppX))
  # Standardize Y
  Y_standardized <- (ppY - y_mean) / y_sd
  
  # Map X to a discrete grid (1:grid_size)
  X_mapped <- round(X_rescaled * (grid_size - 1)) + 1
  
  # Sort by X_mapped
  order_index <- order(X_mapped)
  X_sorted <- X_mapped[order_index]
  Y_sorted <- Y_standardized[order_index]
  
  # Create x_values for plotting (rescale grid to original X range)
  x_values_unit <- seq(0, 1, length.out = grid_size)
  x_values <- x_values_unit * (max(ppX) - min(ppX)) + min(ppX)
  
  return(list(X_sorted = X_sorted,
              Y_sorted = Y_sorted,
              x_values = x_values,
              y_mean = y_mean,
              y_sd = y_sd))
}

# Run wavelet estimation for a given privacy budget (eps)
# eps_values should be a named vector, e.g. c(nonprivate = Inf, eps2 = 2, eps10 = 10)
estimate_wavelet_signals <- function(X_sorted, Y_sorted, grid_size, wavelet_family, boundary, S, s, L_max_np, eps_values) {
  # Compute non-private estimator if needed
  signals <- list()
  
  # For private estimation, compute clipping and sensitivity parameters.
  # We compute the maximum of the mother wavelet at level 0, position 1.
  c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size,
                              family = wavelet_family, bc = boundary, filter_number = S))
  tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max_np)
  sens <- sensitivity(c_psi, length(Y_sorted), tau, log2(grid_size))
  
  for (eps_name in names(eps_values)) {
    eps_val <- eps_values[[eps_name]]
    
    # Choose parameters based on eps
    if (is.infinite(eps_val)) {
      current_tau <- Inf
      current_sens <- 0
    } else {
      current_tau <- tau
      current_sens <- sens
    }
    
    # Estimate detail coefficients
    coeffD_matrix <- private_estimate_coeffD(
      Y = Y_sorted,
      X = X_sorted,
      grid_size = grid_size,
      max_level = log2(grid_size),
      L_max = L_max_np,
      wavelet_family = wavelet_family,
      boundary = boundary,
      S = S,
      tau = current_tau,
      sens = current_sens,
      eps = eps_val
    )
    
    # Estimate coarse coefficient
    coeffC <- private_estimate_coeffC(
      Y = Y_sorted,
      X = X_sorted,
      grid_size = grid_size,
      wavelet_family = wavelet_family,
      boundary = boundary,
      S = S,
      tau = current_tau,
      sens = current_sens,
      eps = eps_val
    )
    
    # Reconstruct the signal on the standardized scale
    est_signal_std <- grid_size * inverse_wavelet(
      coarseC = coeffC,
      coeffD_matrix = coeffD_matrix,
      max_level = log2(grid_size),
      family = wavelet_family,
      filter_number = S,
      bc = boundary
    )
    signals[[eps_name]] <- est_signal_std
  }
  return(signals)
}

# Run the complete analysis pipeline for a given (ppX, ppY)
run_analysis <- function(ppX, ppY, num_bins, grid_size, wavelet_family, boundary, S, s, eps_values) {
  # Uniformly subsample the data (if desired)
  uniform_data <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = num_bins)
  ppX_sub <- uniform_data$ppX
  ppY_sub <- uniform_data$ppY
  
  # Preprocess the data: rescaling, standardization, grid mapping, sorting
  preproc <- preprocess_data(ppX_sub, ppY_sub, grid_size)
  X_sorted <- preproc$X_sorted
  Y_sorted <- preproc$Y_sorted
  x_values <- preproc$x_values
  y_mean <- preproc$y_mean
  y_sd <- preproc$y_sd
  
  # Determine L_max for non-private estimation based on sample size and smoothness parameter s
  L_max_np <- ceiling(log2(length(Y_sorted)) / (2 * s + 1))
  
  # Estimate signals for each eps value
  est_signals_std <- estimate_wavelet_signals(X_sorted, Y_sorted, grid_size,
                                              wavelet_family, boundary, S, s, L_max_np, eps_values)
  # Rescale signals back to the original Y scale
  est_signals <- lapply(est_signals_std, function(sig_std) sig_std * y_sd + y_mean)
  
  return(list(x_values = x_values, 
              ppX_sub = ppX_sub, 
              ppY_sub = ppY_sub,
              est_signals = est_signals))
}

# Plot the original data and estimated signals using ggplot2
plot_wavelet_estimates_ggplot <- function(x_values, ppX, ppY, est_signals, title) {
  # Create a data frame for the original data points
  data_points <- data.frame(x = ppX, y = ppY)
  
  # Create a data frame for the estimated signals
  signal_df <- do.call(rbind, lapply(names(est_signals), function(label) {
    data.frame(x = x_values,
               y = est_signals[[label]],
               privacy = label)
  }))
  
  # Create the ggplot
  p <- ggplot() +
    geom_point(data = data_points, aes(x = x, y = y), color = "black", size = 1) +
    geom_line(data = signal_df, aes(x = x, y = y, color = privacy), size = 1) +
    labs(title = title, x = "Age", y = "Thalach (Original Scale)", color = "Privacy Budget") +
    theme_minimal()
  
  return(p)
}

###############################
# 2. Main Analysis for Each Source
###############################

# Specify parameters
num_bins <- 2        # for uniform subsampling
grid_size <- 512     # grid size (power of 2)
wavelet_family <- "DaubExPhase"
boundary <- "interval"
S <- 3               # filter number
s <- 3               # smoothness parameter
# Define privacy budgets: non-private (Inf), eps = 2, and eps = 10
eps_values <- c(nonprivate = Inf, eps2 = 2, eps10 = 10)

# List of unique sources in the data
sources <- unique(data_filtered$source)

# Create an empty list to store plots for each source
plot_list <- list()

# Run analysis for each source separately and generate ggplot objects
for (src in sources) {
  cat("Running analysis for source:", src, "\n")
  subset_data <- subset(data_filtered, source == src)
  res <- run_analysis(subset_data$age, subset_data$thalach,
                      num_bins, grid_size, wavelet_family, boundary, S, s, eps_values)
  
  p <- plot_wavelet_estimates_ggplot(res$x_values, res$ppX_sub, res$ppY_sub,
                                     res$est_signals, 
                                     title = paste("Wavelet Estimator for Source:", src))
  plot_list[[src]] <- p
}

# Arrange the four plots in a 2x2 grid.
# (If the number of sources is different, adjust the layout accordingly.)
grid.arrange(grobs = plot_list, ncol = 2, nrow = 2)




# -----------------------------
# Additional Code: Train-Test Split and R^2 for Each Source
# -----------------------------

# Define a helper function to compute R^2 using linear interpolation.
compute_R2 <- function(est_signal, x_values, true_x, true_y) {
  # Predict estimated signal at true_x using linear interpolation.
  preds <- approx(x = x_values, y = est_signal, xout = true_x)$y
  ss_res <- sum((true_y - preds)^2, na.rm =T)
  ss_tot <- sum((true_y - mean(true_y))^2, na.rm = T)
  R2 <- 1 - ss_res / ss_tot
  return(R2)
}

# Set seed for reproducibility.
set.seed(123)

# List to store R^2 results for each source.
source_R2_results <- list()

# Loop over each unique source.
for (src in unique(data_filtered$source)) {
  cat("Processing source:", src, "\n")
  
  # Subset the data for the current source.
  source_data <- subset(data_filtered, source == src)
  
  # Train-Test Split (70% train, 30% test) for this source.
  train_idx <- sample(1:nrow(source_data), size = 0.8 * nrow(source_data))
  train_data <- source_data[train_idx, ]
  test_data  <- source_data[-train_idx, ]
  
  # ---- Process Training Data ----
  # Uniformly subsample the training data.
  train_uniform <- subsample_uniform_ppX_ppY(train_data$age, train_data$thalach, num_bins)
  
  # Preprocess the training data.
  train_preproc <- preprocess_data(train_uniform$ppX, train_uniform$ppY, grid_size)
  X_sorted    <- train_preproc$X_sorted
  Y_sorted    <- train_preproc$Y_sorted
  x_values    <- train_preproc$x_values
  y_mean_train <- train_preproc$y_mean
  y_sd_train   <- train_preproc$y_sd
  
  # Determine L_max for non-private estimation based on sample size and smoothness parameter s.
  L_max_np <- ceiling(log2(length(Y_sorted)) / (2 * s + 1))
  
  # Estimate signals for each privacy budget using the training data.
  est_signals_std <- estimate_wavelet_signals(X_sorted, Y_sorted, grid_size,
                                              wavelet_family, boundary, S, s, L_max_np, eps_values)
  # Rescale signals back to original Y scale using training data parameters.
  est_signals <- lapply(est_signals_std, function(sig_std) sig_std * y_sd_train + y_mean_train)
  
  # ---- Compute R^2 on Training and Test Sets ----
  R2_train <- list()
  R2_test  <- list()
  
  for (pb_key in names(est_signals)) {
    # Estimated signal from training.
    est_sig <- est_signals[[pb_key]]
    
    # Compute R^2 on the training set.
    R2_train[[pb_key]] <- compute_R2(est_sig, x_values,
                                     true_x = train_data$age,
                                     true_y = train_data$thalach)
    # Compute R^2 on the test set.
    R2_test[[pb_key]] <- compute_R2(est_sig, x_values,
                                    true_x = test_data$age,
                                    true_y = test_data$thalach)
  }
  
  # Save the results for this source.
  source_R2_results[[src]] <- list(train = R2_train, test = R2_test)
}

# Print the R^2 results for each source.
print(source_R2_results)


# -----------------------------
# Aggregate and Print R^2 Results in a Table
# -----------------------------

# Initialize an empty data frame to store results.
R2_table <- data.frame(Source = character(),
                       Privacy = character(),
                       Train_R2 = numeric(),
                       Test_R2 = numeric(),
                       stringsAsFactors = FALSE)

# Loop over each source in the results.
for (src in names(source_R2_results)) {
  # For each privacy level within the current source.
  for (privacy in names(source_R2_results[[src]]$train)) {
    R2_table <- rbind(R2_table, data.frame(
      Source = src,
      Privacy = privacy,
      Train_R2 = source_R2_results[[src]]$train[[privacy]],
      Test_R2  = source_R2_results[[src]]$test[[privacy]]
    ))
  }
}

# Print the resulting table.
print(R2_table)

