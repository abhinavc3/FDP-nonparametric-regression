library(ggplot2)
library(RColorBrewer)

# Laplace noise generator function
#' Generate Laplace noise.
#'
#' @param n Number of Laplace draws to generate.
#' @param scale Scale parameter of the Laplace distribution.
#' @return A numeric vector of independent Laplace random variables.
rlaplace <- function(n, scale = 1) {
  u <- runif(n, min = -0.5, max = 0.5)
  return(-scale * sign(u) * log(1 - 2 * abs(u)))
}

# Function to compute the private mean with clipping
#' Compute a clipped differentially private mean.
#'
#' @param x Numeric vector to summarize.
#' @param epsilon Privacy budget used for the release.
#' @param tau Symmetric clipping threshold applied before averaging.
#' @return A numeric scalar private mean estimate.
private_mean_clipped <- function(x, epsilon, tau) {
  # Clip the data to the range [-tau, tau]
  x_clipped <- pmin(pmax(x, -tau), tau)
  n <- length(x_clipped)

  # Sensitivity is (2 * tau) / n for the mean (range of clipped values)
  sensitivity <- (2 * tau) / n
  true_mean <- mean(x_clipped)
  noise <- rlaplace(1, sensitivity / epsilon)
  private_mean <- true_mean + noise
  return(private_mean)
}

# Function to compute the private standard deviation with clipping
#' Compute a clipped differentially private standard deviation.
#'
#' @param x Numeric vector to summarize.
#' @param epsilon Privacy budget used for the release.
#' @param tau Symmetric clipping threshold applied before computing the spread.
#' @return A nonnegative numeric private standard deviation estimate.
private_sd_clipped <- function(x, epsilon, tau) {
  # Clip the data to the range [-tau, tau]
  x_clipped <- pmin(pmax(x, -tau), tau)
  n <- length(x_clipped)

  # Sensitivity is (2 * tau) / sqrt(n) for the standard deviation
  sensitivity <- (2 * tau) / sqrt(n)
  true_sd <- sd(x_clipped)
  noise <- rlaplace(1, sensitivity / epsilon)
  private_sd <- max(0, true_sd + noise) # Ensure non-negative result
  return(private_sd)
}

#' Plot private and non-private wavelet fits for a padded one-dimensional dataset.
#'
#' @param ppX Numeric covariate vector.
#' @param ppY Numeric response vector.
#' @param s Smoothness parameter used by the wavelet estimator.
#' @param eps Privacy budget or vector of privacy budgets to plot.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Initial clipping threshold used in private standardization.
#' @param padding_size_percent Fraction of points mirrored at each boundary.
#' @param x_lab_name X-axis label for the plot.
#' @param y_lab_name Y-axis label for the plot.
#' @param grid_size Number of grid points used for wavelet reconstruction.
#' @param wavelet_family Wavelet family passed to the estimator.
#' @param boundary Boundary condition used for wavelet reconstruction.
#' @param include_scatter Logical flag indicating whether raw observations are shown.
#' @param ylim_range Optional y-axis limits passed to `coord_cartesian()`.
#' @param symbol_spacing Step size used when subsampling plotted symbols.
#' @return A `ggplot2` object comparing private and non-private fitted curves.
plot_padded_curve <- function(ppX, ppY, s = 2, eps = 1, S = 4, tau = 6, padding_size_percent = 0.1, x_lab_name = "X", y_lab_name = "Y",
                              grid_size = 2^10, wavelet_family = "DaubExPhase", boundary = "interval",
                              include_scatter = TRUE, ylim_range = NULL, symbol_spacing = 20) {
  # Calculate padding size
  padding_size <- floor(padding_size_percent * length(ppX))

  # Sort ppX and ppY to maintain association
  sorted_indices <- order(ppX)
  ppX_sorted <- ppX[sorted_indices]
  ppY_sorted <- ppY[sorted_indices]

  if (padding_size > 0) {
    # Lower padding (mirror the first 'padding_size' points)
    ppX_lower_padding <- 2 * min(ppX_sorted) - ppX_sorted[padding_size:1]
    ppY_lower_padding <- ppY_sorted[padding_size:1]

    # Upper padding (mirror the last 'padding_size' points)
    ppX_upper_padding <- 2 * max(ppX_sorted) - ppX_sorted[(length(ppX_sorted) - padding_size + 1):length(ppX_sorted)]
    ppY_upper_padding <- ppY_sorted[(length(ppY_sorted) - padding_size + 1):length(ppY_sorted)]

    # Combine padded and original data
    ppX_padded <- c(ppX_lower_padding, ppX_sorted, ppX_upper_padding)
    ppY_padded <- c(ppY_lower_padding, ppY_sorted, ppY_upper_padding)
  } else {
    # No padding applied
    ppX_padded <- ppX_sorted
    ppY_padded <- ppY_sorted
  }

  # Calculate non-private mean and standard deviation
  Y_mean_padded <- mean(ppY_padded)
  Y_sd_padded <- sd(ppY_padded)

  # Standardize ppY using non-private estimates for the non-private estimator
  Y_standardized <- (ppY_padded - Y_mean_padded) / Y_sd_padded

  # Set parameters for wavelet estimation
  ppX_min_padded <- min(ppX_padded)
  ppX_max_padded <- max(ppX_padded)
  ppX_unit <- (ppX_padded - ppX_min_padded) / (ppX_max_padded - ppX_min_padded)
  X <- round(grid_size * ppX_unit)
  print(X)
  max_level <- log2(grid_size)

  # Calculate c_psi
  c_psi <- max(mother_wavelet(
    level = 0, position = 1, grid_size = grid_size,
    family = wavelet_family, bc = boundary, filter_number = S
  ))

  # Ensure no NA/NaN/Inf values in inputs
  if (any(is.na(X)) || any(is.na(Y_standardized))) {
    stop("Error: NA values detected in the inputs for wavelet estimation.")
  }

  # Non-private wavelet estimator
  m <- 1 # Number of servers (if federated)
  ns <- rep(length(Y_standardized), m)
  estimated_signal <- federated_estimator(
    ns, c(Inf), s,
    Y = Y_standardized, X = X, grid_size = grid_size,
    max_level = max_level, wavelet_family = wavelet_family, boundary = boundary,
    S = S, tau = Inf, c_psi = c_psi
  )

  # Reverse the standardization of Y for the non-private estimate
  estimated_signal_original <- estimated_signal * Y_sd_padded + Y_mean_padded

  # Create a data frame for wavelet estimates
  x_values_original <- seq(0, 1, length.out = grid_size) * (ppX_max_padded - ppX_min_padded) + ppX_min_padded

  # Initialize lists to store results
  signals_list <- list()
  line_types_list <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  colors_list <- RColorBrewer::brewer.pal(min(length(eps) + 1, 8), "Dark2")

  # Add the non-private estimate to the list
  signals_list[[1]] <- data.frame(
    x = x_values_original,
    y = estimated_signal_original,
    method = "non-private",
    line_type = "solid",
    color = "black"
  )

  # Iterate over each epsilon value for private estimates
  for (i in seq_along(eps)) {
    current_eps <- eps[i]
    pb <- rep(current_eps, m)
    L_max <- compute_L_max(ns, pb, s)
    tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max)

    # Compute private mean and standard deviation
    Y_private_mean <- private_mean_clipped(ppY_padded, current_eps, tau)
    Y_private_sd <- private_sd_clipped(ppY_padded, current_eps, tau)

    # Standardize ppY using private estimates
    Y_standardized_private <- (ppY_padded - Y_private_mean) / Y_private_sd

    # Private wavelet estimator
    estimated_signal_private <- federated_estimator(
      ns, pb, s,
      Y = Y_standardized_private, X = X, grid_size = grid_size,
      max_level = max_level, wavelet_family = wavelet_family, boundary = boundary,
      S = S, tau = tau, c_psi = c_psi
    )

    # Reverse the standardization of Y for the private estimate
    estimated_signal_private_original <- estimated_signal_private * Y_private_sd + Y_private_mean

    # Add the private estimate to the list
    signals_list[[i + 1]] <- data.frame(
      x = x_values_original,
      y = estimated_signal_private_original,
      method = paste0("eps = ", current_eps),
      line_type = line_types_list[(i %% length(line_types_list)) + 1],
      color = colors_list[(i %% length(colors_list)) + 1]
    )
  }

  # Combine all signals into one data frame
  signals_df <- do.call(rbind, signals_list)

  # Remove estimated values that correspond to the padding
  if (padding_size > 0) {
    unpadded_indices <- (signals_df$x >= min(ppX_sorted)) & (signals_df$x <= max(ppX_sorted))
    signals_df <- signals_df[unpadded_indices, ]
  }

  # Create a data frame for method styles (unique combinations)
  method_styles <- signals_df[!duplicated(signals_df$method), c("method", "color", "line_type")]

  # Create named vectors for scales
  colors <- setNames(method_styles$color, method_styles$method)
  linetypes <- setNames(method_styles$line_type, method_styles$method)
  shapes <- setNames(c(16, 17, 18, 15)[1:length(method_styles$method)], method_styles$method)

  # Subsample for symbols to avoid crowding
  subsampled_df <- signals_df[seq(1, nrow(signals_df), by = symbol_spacing), ]

  # Create the plot
  p <- ggplot()

  # Add scatter plot if requested
  if (include_scatter) {
    p <- p + geom_point(
      data = data.frame(ppX = ppX, ppY = ppY), aes(x = ppX, y = ppY),
      color = "gray", size = 1
    )
  }

  # Map all aesthetics to 'method' in both geoms
  p <- p +
    geom_line(data = signals_df, aes(x = x, y = y, color = method, linetype = method, shape = method), size = 1) +
    geom_point(data = subsampled_df, aes(x = x, y = y, color = method, linetype = method, shape = method), size = 2.5) +
    # Titles and labels
    labs(
      title = "", x = x_lab_name, y = y_lab_name,
      linetype = "Method", color = "Method", shape = "Method"
    ) +
    # Customize scales
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_shape_manual(values = shapes) +
    # Customize theme
    theme_minimal() +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )

  # Apply custom y-axis limits if provided
  if (!is.null(ylim_range)) {
    p <- p + coord_cartesian(ylim = ylim_range)
  }

  return(p)
}


# Function to subsample ppX and ppY for uniform distribution of ppX in [0, 1]
#' Subsample observations so the covariate is represented more uniformly.
#'
#' @param ppX Numeric covariate vector.
#' @param ppY Numeric response vector aligned with `ppX`.
#' @param num_bins Number of bins used to partition the normalized covariate range.
#' @param draw_per_bin Optional cap on the number of sampled observations per bin.
#' @return A list with elements `ppX` and `ppY` containing the subsampled observations.
subsample_uniform_ppX_ppY <- function(ppX, ppY, num_bins, draw_per_bin = NULL) {
  # Normalize ppX to [0, 1] range
  ppX_unit <- (ppX - min(ppX)) / (max(ppX) - min(ppX))

  # Create bins for ppX_unit
  ppX_unit_bins <- cut(ppX_unit, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)

  # Count the number of observations in each bin
  bin_counts <- table(ppX_unit_bins)

  # Determine the maximum number of samples per bin
  max_samples_per_bin <- min(bin_counts)

  # If draw_per_bin is specified, use it; otherwise, use max_samples_per_bin
  if (!is.null(draw_per_bin)) {
    # Use the smaller of draw_per_bin and max_samples_per_bin
    samples_per_bin <- min(draw_per_bin, max_samples_per_bin)
    if (draw_per_bin > max_samples_per_bin) {
      warning("draw_per_bin exceeds the minimum bin count. Using ", samples_per_bin, " samples per bin instead.")
    }
  } else {
    samples_per_bin <- max_samples_per_bin
  }

  # Initialize an empty vector to store the subsample indices
  subsample_indices <- c()

  # Loop over each bin and sample
  for (bin in levels(ppX_unit_bins)) {
    # Get indices of observations in the current bin
    bin_indices <- which(ppX_unit_bins == bin)

    # Randomly sample observations from the bin
    sampled_indices <- sample(bin_indices, samples_per_bin)

    # Append to the subsample indices
    subsample_indices <- c(subsample_indices, sampled_indices)
  }

  # Subsample ppX and ppY using the selected indices
  ppX_subsampled <- ppX[subsample_indices]
  ppY_subsampled <- ppY[subsample_indices]

  return(list(ppX = ppX_subsampled, ppY = ppY_subsampled))
}



#' Plot private and non-private wavelet fits for two datasets on a shared axis.
#'
#' @param ppX1 Numeric covariate vector for the first dataset.
#' @param ppY1 Numeric response vector for the first dataset.
#' @param ppX2 Numeric covariate vector for the second dataset.
#' @param ppY2 Numeric response vector for the second dataset.
#' @param s Smoothness parameter used by the wavelet estimator.
#' @param epsilons Privacy budgets to visualize for each dataset.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Initial clipping threshold used in private standardization.
#' @param padding_size_percent Fraction of points mirrored at each boundary.
#' @param x_lab_name X-axis label for the plot.
#' @param y_lab_name Y-axis label for the plot.
#' @param grid_size Number of grid points used for wavelet reconstruction.
#' @param wavelet_family Wavelet family passed to the estimator.
#' @param boundary Boundary condition used for wavelet reconstruction.
#' @param include_scatter Logical flag indicating whether raw observations are shown.
#' @param ylim_range Optional y-axis limits passed to `coord_cartesian()`.
#' @param symbol_spacing Step size used when subsampling plotted symbols.
#' @param dataset1_name Label used for the first dataset in the legend.
#' @param dataset2_name Label used for the second dataset in the legend.
#' @return A `ggplot2` object comparing fitted curves from the two datasets.
plot_combined_curves_two_datasets <- function(
  ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(0.1, 0.5, 1), S = 4, tau = 6,
  padding_size_percent = 0.1, x_lab_name = "X", y_lab_name = "Y",
  grid_size = 2^10, wavelet_family = "DaubExPhase", boundary = "interval",
  include_scatter = TRUE, ylim_range = NULL, symbol_spacing = 20,
  dataset1_name = "Dataset 1", dataset2_name = "Dataset 2"
) {
  # Define shapes list here to ensure visibility in ggplot components
  shapes_list <- c(16, 17, 18, 15, 4, 8) # Different shapes for the symbols

  # Helper function to process each dataset and epsilon
  process_dataset <- function(ppX, ppY, dataset_name) {
    # Calculate padding size
    padding_size <- floor(padding_size_percent * length(ppX))

    # Sort ppX and ppY to maintain association
    sorted_indices <- order(ppX)
    ppX_sorted <- ppX[sorted_indices]
    ppY_sorted <- ppY[sorted_indices]

    # Create lower and upper padding if needed
    if (padding_size > 0) {
      ppX_lower_padding <- 2 * min(ppX_sorted) - ppX_sorted[padding_size:1]
      ppY_lower_padding <- ppY_sorted[padding_size:1]
      ppX_upper_padding <- 2 * max(ppX_sorted) - ppX_sorted[(length(ppX_sorted) - padding_size + 1):length(ppX_sorted)]
      ppY_upper_padding <- ppY_sorted[(length(ppY_sorted) - padding_size + 1):length(ppY_sorted)]

      # Combine padded and original data
      ppX_padded <- c(ppX_lower_padding, ppX_sorted, ppX_upper_padding)
      ppY_padded <- c(ppY_lower_padding, ppY_sorted, ppY_upper_padding)
    } else {
      ppX_padded <- ppX_sorted
      ppY_padded <- ppY_sorted
    }

    # Calculate non-private mean and standard deviation
    Y_mean_padded <- mean(ppY_padded)
    Y_sd_padded <- sd(ppY_padded)
    Y_standardized <- (ppY_padded - Y_mean_padded) / Y_sd_padded

    # Set parameters for wavelet estimation
    ppX_min_padded <- min(ppX_padded)
    ppX_max_padded <- max(ppX_padded)
    ppX_unit <- (ppX_padded - ppX_min_padded) / (ppX_max_padded - ppX_min_padded)
    X <- round(grid_size * ppX_unit)
    max_level <- log2(grid_size)
    c_psi <- max(mother_wavelet(
      level = 0, position = 1, grid_size = grid_size,
      family = wavelet_family, bc = boundary, filter_number = S
    ))

    # Initialize lists for the results
    signals_list <- list()
    line_types_list <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    colors_list <- RColorBrewer::brewer.pal(min(length(epsilons) + 1, 8), "Dark2")

    # Non-private estimation
    estimated_signal <- federated_estimator(
      rep(length(Y_standardized), 1), c(Inf), s,
      Y = Y_standardized, X = X,
      grid_size = grid_size, max_level = max_level, wavelet_family = wavelet_family,
      boundary = boundary, S = S, tau = Inf, c_psi = c_psi
    ) * Y_sd_padded + Y_mean_padded

    # Generate x values for the estimated signal
    x_values_original <- seq(0, 1, length.out = grid_size) *
      (ppX_max_padded - ppX_min_padded) + ppX_min_padded

    # Add non-private signal to the list
    signals_list[[1]] <- data.frame(
      x = x_values_original,
      y = estimated_signal,
      method = paste0(dataset_name, " non-private"),
      line_type = "solid",
      color = "black",
      shape = shapes_list[1]
    )

    # Private estimates for each epsilon value
    for (i in seq_along(epsilons)) {
      current_eps <- epsilons[i]
      pb <- rep(current_eps, 1)
      L_max <- compute_L_max(rep(length(Y_standardized), 1), pb, s)
      tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max)

      Y_private_mean <- private_mean_clipped(ppY_padded, current_eps, tau)
      Y_private_sd <- private_sd_clipped(ppY_padded, current_eps, tau)
      Y_standardized_private <- (ppY_padded - Y_private_mean) / Y_private_sd

      estimated_signal_private <- federated_estimator(
        rep(length(Y_standardized_private), 1), pb, s,
        Y = Y_standardized_private, X = X,
        grid_size = grid_size, max_level = max_level, wavelet_family = wavelet_family,
        boundary = boundary, S = S, tau = tau, c_psi = c_psi
      ) * Y_private_sd + Y_private_mean

      # Add the private estimate to the list
      signals_list[[i + 1]] <- data.frame(
        x = x_values_original,
        y = estimated_signal_private,
        method = paste0(dataset_name, " eps = ", current_eps),
        line_type = line_types_list[(i %% length(line_types_list)) + 1],
        color = colors_list[(i %% length(colors_list)) + 1],
        shape = shapes_list[(i %% length(shapes_list)) + 1]
      )
    }

    # Combine all signals into one data frame
    signals_df <- do.call(rbind, signals_list)

    # Remove estimated values that correspond to the padding
    if (padding_size > 0) {
      unpadded_indices <- (signals_df$x >= min(ppX_sorted)) & (signals_df$x <= max(ppX_sorted))
      signals_df <- signals_df[unpadded_indices, ]
    }

    return(signals_df)
  }

  # Process both datasets
  signals_df1 <- process_dataset(ppX1, ppY1, dataset1_name)
  signals_df2 <- process_dataset(ppX2, ppY2, dataset2_name)

  # Combine the signals from both datasets
  combined_signals_df <- rbind(signals_df1, signals_df2)

  # Subsample for symbols to avoid crowding
  subsampled_df <- combined_signals_df[seq(1, nrow(combined_signals_df), by = symbol_spacing), ]

  # Create the plot
  p <- ggplot()

  # Add scatter plots if requested with legends for datasets
  if (include_scatter) {
    p <- p +
      geom_point(
        data = data.frame(ppX = ppX1, ppY = ppY1), aes(x = ppX, y = ppY, color = dataset1_name),
        size = 1, shape = 16, alpha = 0.5
      ) +
      geom_point(
        data = data.frame(ppX = ppX2, ppY = ppY2), aes(x = ppX, y = ppY, color = dataset2_name),
        size = 1, shape = 17, alpha = 0.5
      )
  }

  # Add line plots for all estimates
  p <- p +
    geom_line(
      data = combined_signals_df,
      aes(x = x, y = y, linetype = method, color = method), size = 1
    ) +
    geom_point(
      data = subsampled_df,
      aes(x = x, y = y, shape = method, color = method), size = 2.5
    ) +
    scale_shape_manual(values = shapes_list) + # Customize the shapes
    labs(
      title = "", x = x_lab_name, y = y_lab_name, linetype = "Method",
      color = "Dataset/Method", shape = "Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.95, 1.35),
      legend.justification = c("right", "top"),
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14), # Increase x-axis tick label size
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14), # Increase legend text size
      legend.title = element_text(size = 16)
    )

  # Apply custom y-axis limits if provided
  if (!is.null(ylim_range)) {
    p <- p + coord_cartesian(ylim = ylim_range)
  }

  return(p)
}

#' Plot private and non-private wavelet fits for two datasets on a shared axis.
#'
#' @param ppX1 Numeric covariate vector for the first dataset.
#' @param ppY1 Numeric response vector for the first dataset.
#' @param ppX2 Numeric covariate vector for the second dataset.
#' @param ppY2 Numeric response vector for the second dataset.
#' @param s Smoothness parameter used by the wavelet estimator.
#' @param epsilons Privacy budgets to visualize for each dataset.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Initial clipping threshold used in private standardization.
#' @param padding_size_percent Fraction of points mirrored at each boundary.
#' @param x_lab_name X-axis label for the plot.
#' @param y_lab_name Y-axis label for the plot.
#' @param grid_size Number of grid points used for wavelet reconstruction.
#' @param wavelet_family Wavelet family passed to the estimator.
#' @param boundary Boundary condition used for wavelet reconstruction.
#' @param include_scatter Logical flag indicating whether raw observations are shown.
#' @param ylim_range Optional y-axis limits passed to `coord_cartesian()`.
#' @param symbol_spacing Step size used when subsampling plotted symbols.
#' @param dataset1_name Label used for the first dataset in the legend.
#' @param dataset2_name Label used for the second dataset in the legend.
#' @return A `ggplot2` object comparing fitted curves from the two datasets.
plot_combined_curves_two_datasets <- function(
  ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(0.1, 0.5, 1), S = 4, tau = 6,
  padding_size_percent = 0.1, x_lab_name = "X", y_lab_name = "Y",
  grid_size = 2^10, wavelet_family = "DaubExPhase", boundary = "interval",
  include_scatter = TRUE, ylim_range = NULL, symbol_spacing = 20,
  dataset1_name = "Dataset 1", dataset2_name = "Dataset 2"
) {
  # Helper function to process each dataset and epsilon
  process_dataset <- function(ppX, ppY, dataset_name) {
    # Calculate padding size
    padding_size <- floor(padding_size_percent * length(ppX))

    # Sort ppX and ppY to maintain association
    sorted_indices <- order(ppX)
    ppX_sorted <- ppX[sorted_indices]
    ppY_sorted <- ppY[sorted_indices]

    # Create lower and upper padding if needed
    if (padding_size > 0) {
      ppX_lower_padding <- 2 * min(ppX_sorted) - ppX_sorted[padding_size:1]
      ppY_lower_padding <- ppY_sorted[padding_size:1]
      ppX_upper_padding <- 2 * max(ppX_sorted) - ppX_sorted[(length(ppX_sorted) - padding_size + 1):length(ppX_sorted)]
      ppY_upper_padding <- ppY_sorted[(length(ppY_sorted) - padding_size + 1):length(ppY_sorted)]

      # Combine padded and original data
      ppX_padded <- c(ppX_lower_padding, ppX_sorted, ppX_upper_padding)
      ppY_padded <- c(ppY_lower_padding, ppY_sorted, ppY_upper_padding)
    } else {
      ppX_padded <- ppX_sorted
      ppY_padded <- ppY_sorted
    }

    # Calculate non-private mean and standard deviation
    Y_mean_padded <- mean(ppY_padded)
    Y_sd_padded <- sd(ppY_padded)
    Y_standardized <- (ppY_padded - Y_mean_padded) / Y_sd_padded

    # Set parameters for wavelet estimation
    ppX_min_padded <- min(ppX_padded)
    ppX_max_padded <- max(ppX_padded)
    ppX_unit <- (ppX_padded - ppX_min_padded) / (ppX_max_padded - ppX_min_padded)
    X <- round(grid_size * ppX_unit)
    max_level <- log2(grid_size)
    c_psi <- max(mother_wavelet(
      level = 0, position = 1, grid_size = grid_size,
      family = wavelet_family, bc = boundary, filter_number = S
    ))

    # Initialize list for the results
    signals_list <- list()

    # Non-private estimation
    estimated_signal <- federated_estimator(
      rep(length(Y_standardized), 1), c(Inf), s,
      Y = Y_standardized, X = X,
      grid_size = grid_size, max_level = max_level, wavelet_family = wavelet_family,
      boundary = boundary, S = S, tau = Inf, c_psi = c_psi
    ) * Y_sd_padded + Y_mean_padded

    # Generate x values for the estimated signal
    x_values_original <- seq(0, 1, length.out = grid_size) *
      (ppX_max_padded - ppX_min_padded) + ppX_min_padded

    # Add non-private signal to the list
    signals_list[[1]] <- data.frame(
      x = x_values_original,
      y = estimated_signal,
      method = paste0(dataset_name, " non-private"),
      dataset = dataset_name
    )

    # Private estimates for each epsilon value
    for (i in seq_along(epsilons)) {
      current_eps <- epsilons[i]
      pb <- rep(current_eps, 1)
      L_max <- compute_L_max(rep(length(Y_standardized), 1), pb, s)
      tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max)

      Y_private_mean <- private_mean_clipped(ppY_padded, current_eps, tau)
      Y_private_sd <- private_sd_clipped(ppY_padded, current_eps, tau)
      Y_standardized_private <- (ppY_padded - Y_private_mean) / Y_private_sd

      estimated_signal_private <- federated_estimator(
        rep(length(Y_standardized_private), 1), pb, s,
        Y = Y_standardized_private, X = X,
        grid_size = grid_size, max_level = max_level, wavelet_family = wavelet_family,
        boundary = boundary, S = S, tau = tau, c_psi = c_psi
      ) * Y_private_sd + Y_private_mean

      # Add the private estimate to the list
      signals_list[[i + 1]] <- data.frame(
        x = x_values_original,
        y = estimated_signal_private,
        method = paste0(dataset_name, " eps = ", current_eps),
        dataset = dataset_name
      )
    }

    # Combine all signals into one data frame
    signals_df <- do.call(rbind, signals_list)

    # Remove estimated values that correspond to the padding
    if (padding_size > 0) {
      unpadded_indices <- (signals_df$x >= min(ppX_sorted)) & (signals_df$x <= max(ppX_sorted))
      signals_df <- signals_df[unpadded_indices, ]
    }

    return(signals_df)
  }

  # Process both datasets
  signals_df1 <- process_dataset(ppX1, ppY1, dataset1_name)
  signals_df2 <- process_dataset(ppX2, ppY2, dataset2_name)

  # Combine the signals from both datasets
  combined_signals_df <- rbind(signals_df1, signals_df2)

  # Check if combined_signals_df has rows
  if (is.null(nrow(combined_signals_df)) || nrow(combined_signals_df) == 0) {
    stop("Error: combined_signals_df is empty. Check your input data and processing steps.")
  }

  # Ensure symbol_spacing is valid
  max_rows <- nrow(combined_signals_df)
  if (symbol_spacing >= max_rows || symbol_spacing <= 0) {
    symbol_spacing <- max(floor(max_rows / 10), 1) # Adjust to a reasonable value
    warning("symbol_spacing is too large or invalid. Adjusted to ", symbol_spacing)
  }

  # Subsample for symbols to avoid crowding
  subsampled_indices <- seq(1, max_rows, by = symbol_spacing)
  subsampled_df <- combined_signals_df[subsampled_indices, ]

  # Create the plot
  p <- ggplot()

  # Add scatter points first to plot them behind the lines
  if (include_scatter) {
    # Dataset 1 scatter data frame
    scatter_df1 <- data.frame(ppX = ppX1, ppY = ppY1, method = paste0(dataset1_name, " data"))
    # Dataset 2 scatter data frame
    scatter_df2 <- data.frame(ppX = ppX2, ppY = ppY2, method = paste0(dataset2_name, " data"))
    # Combine scatter data frames
    scatter_df <- rbind(scatter_df1, scatter_df2)

    p <- p +
      geom_point(
        data = scatter_df,
        aes(x = ppX, y = ppY, color = method, shape = method, linetype = method),
        size = 0.5, alpha = 0.75
      )
  }

  # Combine method styles from both signals and scatter data
  all_methods <- unique(c(combined_signals_df$method, scatter_df$method))

  # Assign shapes, line types, and colors to methods
  shapes_list <- c(16, 17, 18, 15, 4, 8, 3, 7, 9, 10)
  line_types_list <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  colors_list <- RColorBrewer::brewer.pal(min(length(all_methods), 8), "Dark2")

  # Ensure the lists are long enough
  shapes_list <- rep(shapes_list, length.out = length(all_methods))
  line_types_list <- rep(line_types_list, length.out = length(all_methods))
  colors_list <- rep(colors_list, length.out = length(all_methods))

  # Create a data frame for method styles
  method_styles <- data.frame(
    method = all_methods,
    shape = shapes_list[1:length(all_methods)],
    line_type = ifelse(grepl(" data$", all_methods), "blank", line_types_list[1:length(all_methods)]),
    color = colors_list[1:length(all_methods)],
    stringsAsFactors = FALSE
  )

  # Map styles back to the data frames
  combined_signals_df <- merge(combined_signals_df, method_styles, by = "method")
  subsampled_df <- merge(subsampled_df, method_styles, by = "method")
  if (include_scatter) {
    scatter_df <- merge(scatter_df, method_styles, by = "method")
  }

  # Map all aesthetics to 'method'
  p <- p +
    geom_line(
      data = combined_signals_df,
      aes(x = x, y = y, color = method, linetype = method),
      size = 1
    ) +
    geom_point(
      data = subsampled_df,
      aes(x = x, y = y, color = method, shape = method, linetype = method),
      size = 2.5
    ) +
    # Titles and labels
    labs(
      title = "", x = x_lab_name, y = y_lab_name,
      linetype = "Method", color = "Method", shape = "Method"
    ) +
    # Customize scales
    scale_color_manual(values = setNames(method_styles$color, method_styles$method)) +
    scale_linetype_manual(values = setNames(method_styles$line_type, method_styles$method)) +
    scale_shape_manual(values = setNames(method_styles$shape, method_styles$method)) +
    # Customize theme
    theme_minimal() +
    theme(
      legend.position = c(1, 1), # Move the legend to the top-right
      legend.justification = c("right", "top"), # Align it to the top-right
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )

  # Use guides to adjust legend appearance for scatter points
  p <- p + guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1.0)
    ),
    shape = guide_legend(
      override.aes = list(size = 3, alpha = 1.0)
    )
  )

  # Apply custom y-axis limits if provided
  if (!is.null(ylim_range)) {
    p <- p + coord_cartesian(ylim = ylim_range)
  }

  return(p)
}
