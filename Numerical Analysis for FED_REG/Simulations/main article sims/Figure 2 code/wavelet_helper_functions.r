library(wavethresh)

# clear the workspace
examples <- FALSE

# Function to generate a mother wavelet basis function for a given level and position
#' Construct a mother wavelet basis function on a discrete grid.
#'
#' @param level Detail level in the wavelet decomposition.
#' @param position Position index of the basis function within the level.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param bc Boundary condition passed to `wavethresh::wd()`.
#' @param family Wavelet family passed to `wavethresh::wd()`.
#' @param filter_number Filter number used by the selected wavelet family.
#' @return A numeric vector containing the reconstructed mother wavelet basis function.
mother_wavelet <- function(level, position, grid_size, bc = "interval", family = "DaubExPhase", filter_number = 1) {
  # Access the detail coefficients at the given level
  decomp <- wd(rep(0, grid_size), filter.number = filter_number, family = family, bc = bc)

  coeffs_length <- length(accessD(decomp, level))

  # Ensure the position is within the valid range
  if (position < 1 || position > coeffs_length) {
    stop(paste("Position must be between 1 and", coeffs_length, "for level", level))
  }

  # Set the detail coefficient at the specified position
  coeffs <- rep(0, coeffs_length)
  coeffs[position] <- 1
  decomp <- putD(decomp, level, coeffs)

  # Reconstruct the wavelet basis function
  wavelet_function <- wr(decomp)

  return(wavelet_function)
}

#' Construct a father wavelet basis function on a discrete grid.
#'
#' @param level Approximation level in the wavelet decomposition.
#' @param position Position index of the scaling basis function within the level.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param bc Boundary condition passed to `wavethresh::wd()`.
#' @param family Wavelet family passed to `wavethresh::wd()`.
#' @param filter_number Filter number used by the selected wavelet family.
#' @return A numeric vector containing the reconstructed father wavelet basis function.
father_wavelet <- function(level = 0, position = 1, grid_size, bc = "interval", family = "DaubExPhase", filter_number = 1) {
  if (bc == "periodic") {
    decomp <- wd(c(1, rep(0, grid_size - 1)), filter.number = filter_number, family = family, bc = bc)
  }
  if (bc == "symmetric") {
    decomp <- wd(c(1, rep(0, grid_size - 1)), filter.number = filter_number, family = family, bc = bc)
  }
  if (bc == "interval") {
    # Access the approximation coefficients at the given level
    decomp <- wd(rep(0, grid_size), filter.number = filter_number, family = family, bc = bc)
  } else {
    print("Boundary condition not recognized")
  }

  coeffs_length <- length(accessC(decomp, level))

  # Ensure the position is within the valid range
  if (position < 1 || position > coeffs_length) {
    stop(paste("Position must be between 1 and", coeffs_length, "for level", level))
  }

  # Set the approximation (scaling) coefficient at the specified position
  coeffs <- rep(0, coeffs_length)
  coeffs[position] <- 1
  decomp <- putC(decomp, level, coeffs)

  # Reconstruct the scaling (father wavelet) basis function
  wavelet_function <- wr(decomp)


  return(wavelet_function)
}


# PLOTTING A WAVELET BASIS FUNCTION EXAMPLE
if (examples) {
  # Number of points in the signal (should be a power of 2)
  grid_size <- 2^15
  max_level <- log2(grid_size)
  boundary <- "interval"
  wavelet_family <- "Coiflets"
  # Other wavelet families considered for these examples include DaubLeAsymm, Coiflets, Symmlets, and Haar.
  S <- 3 # 1 for Haar
  x_values <- seq(0, 1, length.out = grid_size)

  l <- 1
  k <- 1
  # Example: Visualize the wavelet basis function at level 3, position 5
  wavelet <- mother_wavelet(level = l, position = k, grid_size = grid_size, family = wavelet_family, bc = boundary, filter_number = S)
  wavelet2 <- mother_wavelet(level = l, position = k + 1, grid_size = grid_size, family = wavelet_family, bc = boundary, filter_number = S)

  l <- 0
  k <- 1
  fwavelet <- father_wavelet(level = l, position = k, grid_size = grid_size, family = wavelet_family, bc = boundary, filter_number = S)

  # Check orthonormality and vanishing 0th moment
  sum(wavelet)
  sum(wavelet^2)
  sum(wavelet * wavelet2)

  quartz(width = 10, height = 7)
  x_values <- seq(0, 1, length.out = grid_size)
  # Plot the wavelet basis function
  plot(
    x = x_values, y = wavelet, type = "l", col = "blue", lwd = 2,
    main = "Daubechies Wavelet Basis Function",
    xlab = "Time", ylab = "Amplitude"
  )
  # plot wavelet2 in red
  lines(x = x_values, y = wavelet2, type = "l", col = "red", lwd = 2)
  lines(x = x_values, y = fwavelet, type = "l", col = "green", lwd = 2)


  coeffD_matrix <- matrix(NA, nrow = grid_size, ncol = max_level)
  for (i in 1:(max_level)) {
    coeffD_matrix[1:2^(i - 1), i] <- 0
  }
}

# Function generating a function from a vector of coefficients at resolution level
#' Reconstruct a signal from coarse and detail wavelet coefficients.
#'
#' @param coarseC Coarse scaling coefficient or vector of scaling coefficients.
#' @param coeffD_matrix Matrix of detail coefficients indexed by level and position.
#' @param max_level Maximum decomposition level to reconstruct.
#' @param bc Boundary condition passed to `wavethresh::wd()`.
#' @param family Wavelet family passed to `wavethresh::wd()`.
#' @param filter_number Filter number used by the selected wavelet family.
#' @return A numeric vector containing the reconstructed signal on the working grid.
inverse_wavelet <- function(coarseC = 0, coeffD_matrix, max_level, bc = "interval", family = "DaubExPhase", filter_number = 1) {
  gs <- dim(coeffD_matrix)[1]
  # Initialize the signal with the coarse coefficients
  decomp <- wd(rep(0, gs), filter.number = filter_number, family = family, bc = bc)

  # Set the coarse coefficients
  decomp <- putC(decomp, 0, coarseC)

  # Set the detail coefficients
  for (i in 1:(max_level)) {
    decomp <- putD(decomp, i - 1, coeffD_matrix[1:2^(i - 1), i])
  }

  # Reconstruct the signal
  signal <- wr(decomp)
  return(signal)
}

if (examples) {
  coeffD_matrix[1, 1] <- 1
  coeffD_matrix[2, 3] <- 10
  wavelet <- inverse_wavelet(coarseC = 0, coeffD_matrix = coeffD_matrix, max_level = max_level, family = wavelet_family, bc = boundary, filter_number = S)

  plot(
    x = x_values, y = wavelet, type = "l", col = "blue", lwd = 2,
    main = "Daubechies Wavelet Basis Function (Level 3, Position 5)",
    xlab = "Time", ylab = "Amplitude"
  )
}

# Function to generate Rademacher variables
#' Sample Rademacher random variables.
#'
#' @param size Number of draws to generate.
#' @return A numeric vector of length `size` with entries in `{-1, 1}`.
rrad <- function(size = 1) {
  return(sample(c(-1, 1), size, replace = TRUE))
}

# Calculate the Besov norm based on the coefficients
#' Compute the Besov norm implied by a wavelet coefficient matrix.
#'
#' @param coeffD_matrix Matrix of wavelet detail coefficients.
#' @param s Smoothness parameter in the Besov space definition.
#' @param p Integrability parameter controlling within-level aggregation.
#' @param q Integrability parameter controlling across-level aggregation.
#' @return A numeric scalar giving the Besov norm of the supplied coefficients.
besov_norm <- function(coeffD_matrix, s, p, q) {
  max_level <- ncol(coeffD_matrix)
  norm_sum <- 0

  for (j in 1:max_level) {
    level_coeffs <- coeffD_matrix[1:2^(j - 1), j]
    inner_sum <- sum(abs(level_coeffs)^p)^(1 / p)
    scaling_factor <- 2^(j * (s + 1 / 2 - 1 / p))

    if (q < Inf) {
      norm_sum <- norm_sum + (scaling_factor * inner_sum)^q
    } else {
      norm_sum <- max(norm_sum, scaling_factor * inner_sum)
    }
  }

  if (q < Inf) {
    return(norm_sum^(1 / q))
  } else {
    return(norm_sum)
  }
}

# Generate a Besov-smooth function
#' Generate a Besov-smooth test signal from random wavelet coefficients.
#'
#' @param grid_size Number of grid points used for the reconstructed signal.
#' @param s Smoothness parameter in the Besov space definition.
#' @param p Integrability parameter controlling within-level aggregation.
#' @param q Integrability parameter controlling across-level aggregation.
#' @param bc Boundary condition used when reconstructing the signal.
#' @return A numeric vector containing the reconstructed Besov-smooth signal.
besov_smooth_function <- function(grid_size, s = 0.5, p = 2, q = 2, bc = "periodic") {
  max_level <- log2(grid_size)
  coeffD_matrix <- matrix(0, nrow = grid_size, ncol = max_level)

  # Populate the coefficient matrix with Rademacher variables
  for (level in 1:max_level) {
    num_coeffs <- 2^(level - 1)
    rademacher_vars <- rrad(num_coeffs)

    # Scaling for Besov smoothness (this can be adjusted based on the desired smoothness)
    scaling_factor <- 2^(-level * (s + 1 / 2)) * 2^(level / p) * (2)^(-level / q)

    # Populate the coefficients with scaled Rademacher variables
    coeffD_matrix[1:num_coeffs, level] <- rademacher_vars * scaling_factor
  }
  print(besov_norm(coeffD_matrix, s, p, q))

  # Reconstruct the signal using the inverse wavelet transform
  reconstructed_signal <- inverse_wavelet(coarseC = 0, coeffD_matrix = coeffD_matrix, max_level = max_level, family = "DaubExPhase", bc = bc, filter_number = S)
}
if (examples) {
  bf <- besov_smooth_function(grid_size = grid_size, s = 1.5, p = 2, q = 2, bc = "interval")
  plot(
    x = x_values, y = bf, type = "l", col = "blue", lwd = 2,
    main = "Besov-Smooth Function",
    xlab = "x (Time)", ylab = "Amplitude"
  )
}

# Define the clipping function
#' Clip values to a symmetric threshold.
#'
#' @param x Numeric vector to be clipped.
#' @param tau Nonnegative clipping threshold; `Inf` leaves the input unchanged.
#' @return A numeric vector with each entry truncated to `[-tau, tau]`.
clip <- function(x, tau = Inf) {
  return(pmin(pmax(x, -tau), tau))
}

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

# Function to estimate the wavelet coefficient \hat{f}^{(j)}_{lk;\tau}
#' Estimate a wavelet coefficient from observations on the working grid.
#'
#' @param Y Numeric response vector.
#' @param X Integer grid locations associated with `Y`.
#' @param psi Numeric vector containing the evaluated wavelet basis function.
#' @param tau Clipping threshold applied to `Y` before estimation.
#' @return A numeric scalar estimate of the requested wavelet coefficient.
estimate_wavelet_coefficient <- function(Y, X, psi, tau = Inf) {
  # Apply the clipping function to Y
  Y_clipped <- clip(Y, tau)

  # Compute the product of the thresholded Y and the wavelet function psi evaluated at X
  products <- Y_clipped * psi[X]
  # Estimate the wavelet coefficient by taking the mean of the products
  estimated_coefficient <- mean(products)

  return(estimated_coefficient)
}

# Calculate the sensitivity of the wavelet coefficient estimation
#' Compute the sensitivity bound used in private coefficient estimation.
#'
#' @param c_psi Upper bound on the absolute value of the mother wavelet.
#' @param nj Number of observations held by one machine or sample block.
#' @param tau Clipping threshold applied to the responses.
#' @param L_max Maximum wavelet resolution level retained in the estimator.
#' @return A numeric scalar sensitivity bound for the corresponding private release.
sensitivity <- function(c_psi, nj, tau, L_max) {
  return(2 * c_psi * tau * sqrt(sum(2^(0:L_max))) / nj)
}

# Function to privately estimate the wavelet coefficients (detail)
#' Privately estimate detail wavelet coefficients for one sample block.
#'
#' @param Y Numeric response vector.
#' @param X Integer grid locations associated with `Y`.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param max_level Maximum decomposition level available on the grid.
#' @param L_max Maximum level actively estimated under the privacy budget.
#' @param wavelet_family Wavelet family used in basis construction.
#' @param boundary Boundary condition used in basis construction.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Clipping threshold applied to the responses.
#' @param sens Sensitivity bound used to scale Laplace noise.
#' @param eps Privacy budget used for the coefficient release.
#' @return A matrix of estimated detail coefficients indexed by level and position.
private_estimate_coeffD <- function(Y, X, grid_size, max_level, L_max, wavelet_family, boundary, S, tau = Inf, sens = 0, eps = Inf) {
  # Initialize the coefficient matrix
  coeffD_matrix <- matrix(0, nrow = grid_size, ncol = max_level)

  # Loop through levels and positions to calculate coefficients
  for (level in 1:L_max) {
    for (position in 1:2^(level - 1)) {
      # Generate the mother wavelet for the current level and position
      psi <- mother_wavelet(level - 1, position, grid_size, family = wavelet_family, bc = boundary, filter_number = S)

      if (eps == Inf) {
        # Estimate the wavelet coefficient using the thresholded values
        estimated_coefficient <- estimate_wavelet_coefficient(Y, X, psi, Inf)
      }

      # Add Laplace noise to ensure privacy
      if (eps < Inf) {
        # Estimate the wavelet coefficient using the thresholded values
        estimated_coefficient <- estimate_wavelet_coefficient(Y, X, psi, tau)

        # Scale the noise by sensitivity and privacy budget (epsilon)
        noise <- rlaplace(1, scale = sens / eps)
        estimated_coefficient <- estimated_coefficient + noise
      }

      # Estimate the wavelet coefficient using the thresholded values
      coeffD_matrix[position, level] <- estimated_coefficient
    }
  }

  return(coeffD_matrix)
}

# Function to privately estimate mean (coarse) coefficient
#' Privately estimate the coarse scaling coefficient for one sample block.
#'
#' @param Y Numeric response vector.
#' @param X Integer grid locations associated with `Y`.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param wavelet_family Wavelet family used in basis construction.
#' @param boundary Boundary condition used in basis construction.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Clipping threshold applied to the responses.
#' @param sens Sensitivity bound used to scale Laplace noise.
#' @param eps Privacy budget used for the coefficient release.
#' @return A numeric scalar estimate of the coarse scaling coefficient.
private_estimate_coeffC <- function(Y, X, grid_size, wavelet_family, boundary, S, tau = Inf, sens = 0, eps = Inf) {
  # Generate the father wavelet function
  phi <- father_wavelet(level = 0, position = 1, grid_size = grid_size, family = wavelet_family, bc = boundary, filter_number = S)

  # Add Laplace noise to ensure privacy
  if (eps < Inf) {
    # Estimate the mean coefficient using the thresholded values
    estimated_coefficient <- estimate_wavelet_coefficient(Y, X, phi, tau)

    # Scale the noise by sensitivity and privacy budget (epsilon)
    noise <- rlaplace(1, scale = sens / eps)
    estimated_coefficient <- estimated_coefficient + noise
  }
  if (eps == Inf) {
    # Estimate the mean coefficient using the thresholded values
    estimated_coefficient <- estimate_wavelet_coefficient(Y, X, phi, Inf)
  }

  return(estimated_coefficient)
}


# Example usage
if (examples) {
  set.seed(123)
  nj <- 100 # Number of observations
  s <- 2
  eps <- 10
  # Private L_max
  L_max <- min(ceiling(log2(nj^2 * eps^2) / (2 * s + 2)), ceiling(log2(nj) / (2 * s + 1)))
  # Non-private L_max
  L_max_np <- ceiling(log2(nj) / (2 * s + 1))

  # Generate uniform draws from 1 to grid_size
  X <- sample(1:grid_size, nj, replace = TRUE)

  # Generate the Besov-smooth function
  signal <- 10 * sqrt(grid_size) * besov_smooth_function(grid_size = grid_size, s = s, p = 2, q = 2, bc = boundary)

  # Generate the noisy observations
  Y <- signal[X] + rnorm(nj, mean = 0, sd = 1)

  # Determine the range of Y-axis
  y_min <- min(Y, signal)
  y_max <- max(Y, signal)

  # Extend the Y-axis range by adding some padding
  y_padding <- 0.2 * (y_max - y_min)
  ylim_range <- c(y_min - y_padding, y_max + y_padding)

  # plot the noisy observations and the true signal
  plot(
    x = x_values, y = signal, type = "l", col = "blue", lwd = 2,
    main = "True Signal vs. Noisy Observations",
    xlab = "x (Time)", ylab = "Amplitude", ylim = ylim_range
  )
  points(x = X / grid_size, y = Y, col = "black", bg = NA, pch = 1)

  # Estimate the wavelet coefficients (detail) non-privately
  coeffD_matrix <- private_estimate_coeffD(Y, X, grid_size, max_level, L_max_np, wavelet_family, boundary, S, tau = Inf, sens = 0, eps = Inf)

  # Estimate the mean coefficient (coarse) non-privately
  coeffC <- private_estimate_coeffC(Y, X, grid_size, wavelet_family, boundary, S, tau = Inf, sens = 0, eps = Inf)

  # Reconstruct the estimated signal using the inverse wavelet transform
  est_signal_non_private <- (grid_size) * inverse_wavelet(coarseC = coeffC, coeffD_matrix = coeffD_matrix, max_level = max_level, family = wavelet_family, filter_number = S, bc = boundary)


  # Calculate the maximum of the mother wavelet
  c_psi <- max(mother_wavelet(level = 0, position = 1, grid_size = grid_size, family = wavelet_family, bc = boundary, filter_number = S))

  # Set the clipping parameter tau
  tau <- sqrt(grid_size) * c_psi + sqrt((2 * s + 1) * L_max_np)

  # Calculate the sensitivity
  sens <- sensitivity(c_psi, nj, tau = tau, L_max = L_max)

  # Private estimation of the wavelet coefficients (detail)
  coeffD_matrix_private <- private_estimate_coeffD(Y, X, grid_size, max_level, L_max, wavelet_family, boundary, S, tau = tau, sens = sens, eps = eps)

  # Private estimation of the mean coefficient (coarse)
  coeffC_private <- private_estimate_coeffC(Y, X, grid_size, wavelet_family, boundary, S, tau = tau, sens = sens, eps = eps)

  # Reconstruct the estimated signal using the inverse wavelet transform
  est_signal_private <- (grid_size) * inverse_wavelet(coarseC = coeffC_private, coeffD_matrix = coeffD_matrix_private, max_level = max_level, family = wavelet_family, filter_number = S, bc = boundary)

  # Plot the original and reconstructed signals, together with the observations
  plot(
    x = x_values, y = signal, type = "l", col = "blue", lwd = 2,
    main = "Original vs. Reconstructed Signal",
    xlab = "x (Time)", ylab = "Amplitude", ylim = ylim_range
  )
  lines(x = x_values, y = est_signal_non_private, col = "red", lwd = 2, lty = 2)
  lines(x = x_values, y = est_signal_private, col = "green", lwd = 2, lty = 2)
  points(x = X / grid_size, y = Y, col = "black", bg = NA, pch = 1)
}

## Distributed setup

# Function to compute D based on the equation
#' Solve the fixed-point equation that defines the effective resolution parameter.
#'
#' @param ns Integer vector of local sample sizes across machines.
#' @param pb Numeric vector of privacy budgets across machines.
#' @param s Smoothness parameter in the theoretical rate expressions.
#' @return A numeric scalar solution for the effective resolution parameter `D`.
solve_for_D <- function(ns, pb, s) {
  # Define the function for D based on the given equation
  equation_for_D <- function(D) {
    lhs <- D^(2 * s + 2)
    rhs <- sum(pmin(ns^2 * pb^2, ns * D))
    return(lhs - rhs)
  }

  # Solve for D using uniroot
  D_solution <- uniroot(equation_for_D, interval = c(1, 1e5))$root
  return(D_solution)
}

# Compute the resolution level given the observation vector
#' Compute the maximum wavelet resolution level used by the estimator.
#'
#' @param ns Integer vector of local sample sizes across machines.
#' @param pb Numeric vector of privacy budgets across machines.
#' @param s Smoothness parameter in the theoretical rate expressions.
#' @return An integer giving the largest resolution level retained by the procedure.
compute_L_max <- function(ns, pb, s) {
  if (all(pb == Inf)) {
    # When epsilon is infinite, set L_max based on total observations
    N <- sum(ns)
    L_max <- ceiling(log2(N) / (2 * s + 1))
  } else {
    D <- solve_for_D(ns, pb, s)
    L_max <- max(1, ceiling(log2(D)))
  }
  return(L_max)
}

# Function to compute the weights u_j based on v_j
#' Compute aggregation weights for the federated estimator.
#'
#' @param ns Integer vector of local sample sizes across machines.
#' @param pb Numeric vector of privacy budgets across machines.
#' @param L Resolution level used in the weighting formula.
#' @return A numeric vector of nonnegative aggregation weights summing to one.
compute_weights <- function(ns, pb, L) {
  if (all(pb == Inf)) {
    u_js <- ns / sum(ns)
  } else {
    v_js <- pmin(ns^2 * pb^2, ns * 2^L)
    u_js <- v_js / sum(v_js)
  }
  return(u_js)
}

# Compute federated private estimator for coarse coefficients
#' Aggregate private coarse coefficient estimates across machines.
#'
#' @param ns Integer vector of local sample sizes across machines.
#' @param pb Numeric vector of privacy budgets across machines.
#' @param s Smoothness parameter in the theoretical rate expressions.
#' @param Y Numeric response vector stacked across machines.
#' @param X Integer grid locations stacked across machines.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param max_level Maximum decomposition level available on the grid.
#' @param wavelet_family Wavelet family used in basis construction.
#' @param boundary Boundary condition used in basis construction.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Clipping threshold applied to the responses.
#' @param c_psi Upper bound on the absolute value of the mother wavelet.
#' @return A numeric scalar giving the aggregated coarse scaling coefficient.
federated_private_estimate_coeffC <- function(ns, pb, s, Y, X, grid_size, max_level, wavelet_family, boundary, S, tau = Inf, c_psi = 0) {
  m <- length(ns) # Number of servers (machines)
  L_max <- compute_L_max(ns, pb, s) # Compute the maximum resolution level
  coeffCs <- numeric(m) # Store coarse coefficients from each machine

  start_index <- 1
  # Loop over each machine
  for (j in 1:m) {
    # Get indices for local data for machine j
    end_index <- start_index + ns[j] - 1
    sens <- sensitivity(c_psi, ns[j], tau = tau, L_max = L_max)

    # Local data for machine j
    Yj <- Y[start_index:end_index]
    Xj <- X[start_index:end_index]

    # Update start index for the next machine
    start_index <- end_index + 1

    # Estimate the coarse coefficient privately for machine j
    coeffCs[j] <- private_estimate_coeffC(Yj, Xj, grid_size, wavelet_family, boundary, S, tau = tau, sens = sens, eps = pb[j])
  }

  # Compute the weights for aggregation
  u_js <- compute_weights(ns, pb, L_max)

  # Return the federated coarse coefficient as the average of local estimates
  return(sum(u_js * coeffCs))
}

# Compute federated private estimator for detail coefficients
#' Aggregate private detail coefficient estimates across machines.
#'
#' @param ns Integer vector of local sample sizes across machines.
#' @param pb Numeric vector of privacy budgets across machines.
#' @param s Smoothness parameter in the theoretical rate expressions.
#' @param Y Numeric response vector stacked across machines.
#' @param X Integer grid locations stacked across machines.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param max_level Maximum decomposition level available on the grid.
#' @param wavelet_family Wavelet family used in basis construction.
#' @param boundary Boundary condition used in basis construction.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Clipping threshold applied to the responses.
#' @param c_psi Upper bound on the absolute value of the mother wavelet.
#' @return A matrix of aggregated detail coefficients indexed by level and position.
federated_private_estimate_coeffD <- function(ns, pb, s, Y, X, grid_size, max_level, wavelet_family, boundary, S, tau = Inf, c_psi = 0) {
  m <- length(ns) # Number of servers (machines)
  L_max <- compute_L_max(ns, pb, s) # Compute the maximum resolution level
  coeffDs_list <- vector("list", m) # List to store detail coefficients for each machine

  # Initialize coefficient matrices for each machine
  for (j in 1:m) {
    coeffDs_list[[j]] <- matrix(0, nrow = grid_size, ncol = max_level)
  }

  start_index <- 1
  # Loop over each machine
  for (j in 1:m) {
    # Get indices for local data for machine j
    end_index <- start_index + ns[j] - 1
    sens <- sensitivity(c_psi, ns[j], tau = tau, L_max = L_max)

    # Local data for machine j
    Yj <- Y[start_index:end_index]
    Xj <- X[start_index:end_index]

    # Update start index for the next machine
    start_index <- end_index + 1

    # Estimate the detail coefficients privately for machine j
    coeffDs_list[[j]] <- private_estimate_coeffD(Yj, Xj, grid_size, max_level, L_max = L_max, wavelet_family, boundary, S, tau = tau, sens = sens, eps = pb[j])
  }

  # Compute the weights for aggregation
  u_js <- compute_weights(ns, pb, L_max)

  # Aggregate the detail coefficients across machines
  sum_coefs <- matrix(0, nrow = grid_size, ncol = max_level)
  for (j in 1:m) {
    sum_coefs <- sum_coefs + u_js[j] * coeffDs_list[[j]]
  }

  # Return the federated detail coefficient as the average of local estimates
  return(sum_coefs)
}

# Final federated estimator function
#' Reconstruct the final federated nonparametric estimator.
#'
#' @param ns Integer vector of local sample sizes across machines.
#' @param pb Numeric vector of privacy budgets across machines.
#' @param s Smoothness parameter in the theoretical rate expressions.
#' @param Y Numeric response vector stacked across machines.
#' @param X Integer grid locations stacked across machines.
#' @param grid_size Number of grid points used for the discrete signal.
#' @param max_level Maximum decomposition level available on the grid.
#' @param wavelet_family Wavelet family used in basis construction.
#' @param boundary Boundary condition used in basis construction.
#' @param S Filter number used by the selected wavelet family.
#' @param tau Clipping threshold applied to the responses.
#' @param c_psi Upper bound on the absolute value of the mother wavelet.
#' @return A numeric vector containing the reconstructed federated estimator on the working grid.
federated_estimator <- function(ns, pb, s, Y, X, grid_size, max_level, wavelet_family, boundary, S, tau = Inf, c_psi = 0) {
  coeffC <- 0
  # Estimate the coarse coefficient, privately and federated
  coeffC <- federated_private_estimate_coeffC(ns, pb, s, Y = Y, X = X, grid_size, max_level, wavelet_family, boundary, S, tau = tau, c_psi = c_psi)

  coeffD_matrix <- matrix(0, nrow = grid_size, ncol = max_level)
  # Estimate the detail coefficients, privately and federated
  coeffD_matrix <- federated_private_estimate_coeffD(ns, pb, s, Y = Y, X = X, grid_size, max_level, wavelet_family, boundary, S, tau = tau, c_psi = c_psi)

  # Reconstruct the final signal using the inverse wavelet transform
  estimated_signal <- grid_size * inverse_wavelet(coarseC = coeffC, coeffD_matrix = coeffD_matrix, max_level = max_level, family = wavelet_family, filter_number = S, bc = boundary)

  return(estimated_signal)
}

# Function to calculate RMSE
#' Compute the root mean squared error between two signals.
#'
#' @param original_signal Numeric vector containing the reference signal.
#' @param reconstructed_signal Numeric vector containing the estimated signal.
#' @return A numeric scalar root mean squared error.
calculate_rmse <- function(original_signal, reconstructed_signal) {
  return(sqrt(mean((original_signal - reconstructed_signal)^2)))
}

# Function to calculate MAE
#' Compute the mean absolute error between two signals.
#'
#' @param original_signal Numeric vector containing the reference signal.
#' @param reconstructed_signal Numeric vector containing the estimated signal.
#' @return A numeric scalar mean absolute error.
calculate_mae <- function(original_signal, reconstructed_signal) {
  return(mean(abs(original_signal - reconstructed_signal)))
}

# Function to calculate RMSE at a given point t
#' Compute the pointwise squared error at a target location.
#'
#' @param original_signal Numeric vector containing the reference signal.
#' @param reconstructed_signal Numeric vector containing the estimated signal.
#' @param t Target location at which the pointwise error is evaluated.
#' @param x_values Grid locations associated with the signal vectors.
#' @return A numeric scalar squared error at the grid point closest to `t`.
calculate_mse_at_t <- function(original_signal, reconstructed_signal, t, x_values) {
  # Find the index of the closest point in x_values to the specified t
  index <- which.min(abs(x_values - t))
  # Compute RMSE at this single point
  return((original_signal[index] - reconstructed_signal[index])^2)
}

# Function to calculate MAE at a given point t
#' Compute the pointwise absolute error at a target location.
#'
#' @param original_signal Numeric vector containing the reference signal.
#' @param reconstructed_signal Numeric vector containing the estimated signal.
#' @param t Target location at which the pointwise error is evaluated.
#' @param x_values Grid locations associated with the signal vectors.
#' @return A numeric scalar absolute error at the grid point closest to `t`.
calculate_mae_at_t <- function(original_signal, reconstructed_signal, t, x_values) {
  # Find the index of the closest point in x_values to the specified t
  index <- which.min(abs(x_values - t))
  # Compute MAE at this single point
  return(abs(original_signal[index] - reconstructed_signal[index]))
}

if (examples) {
  # Example usage of RMSE and MAE at a specific point t
  t <- 0.5 # Point at which you want to evaluate the error


  # Example use of the above functions
  m <- 5 # Set the number of servers
  # Set the number of local observations for each machine
  ns <- rep(100, m) # initially set to 1 observation per machine
  N <- sum(ns)
  L_max <- compute_L_max(ns, pb, s)

  # Generate data according to ns
  X <- sample(1:grid_size, N, replace = TRUE)
  Y <- signal[X] + rnorm(N, mean = 0, sd = 0.1)


  # Set the privacy budget for each machine
  pb <- rep(10, m) # initially set to 1

  compute_L_max(ns, pb, s)

  # Assuming ns, pb, Y, X, grid_size, etc. are already defined
  estimated_signal_private <- federated_estimator(ns, pb, s, Y = Y, X = X, grid_size, max_level, wavelet_family, boundary, S, tau = tau, c_psi = c_psi)

  estimated_signal_non_private <- federated_estimator(rep(N, 1), rep(Inf, 1), s, Y = Y, X = X, grid_size, max_level, wavelet_family, boundary, S, tau = Inf, c_psi = 1)

  # Plot the original and reconstructed signals, together with the observations
  plot(
    x = x_values, y = signal, type = "l", col = "blue", lwd = 2,
    main = "Original vs. Reconstructed Signal",
    xlab = "x (Time)", ylab = "Amplitude", ylim = ylim_range
  )
  lines(x = x_values, y = estimated_signal_non_private, col = "red", lwd = 2, lty = 2)
  lines(x = x_values, y = estimated_signal_private, col = "green", lwd = 2, lty = 2)
  points(x = X / grid_size, y = Y, col = "black", bg = NA, pch = 1)
}
