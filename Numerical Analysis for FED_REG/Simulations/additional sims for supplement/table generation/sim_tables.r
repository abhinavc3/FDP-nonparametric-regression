# Load necessary libraries
library(wavethresh)
if (!require(xtable)) install.packages("xtable")
library(xtable)
# Set working directory
#setwd("/Users/lassev/Dropbox (Penn)/Distributed comm. and privacy constraints/Federated Learning for Nonparametric Regression/Simulations/")
setwd("/home/stat/lassev/PrivateNonparametricRegression")

# Include wavelet helper functions
source("wavelet_helper_functions.r")

# Set seed for reproducibility
set.seed(2024)

# Simulation parameters
epsilon_values <- c(0.1,0.5,1,2,5,10,Inf)
machine_counts <- c(1,100,250,2000)
N <- 2000  # Total number of observations
num_simulations <- 1000  # Number of simulations per configuration
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

# Initialize results lists
results_mse_total <- list()
results_mse_point <- list()

# Function to calculate index of point t
get_index_t <- function(t_point, x_values) {
  which.min(abs(x_values - t_point))
}

index_t <- get_index_t(t_point, x_values)

# Loop over each epsilon value
for (eps in epsilon_values) {
  # Print progress
  cat("Epsilon:", eps, "\n")

  # Initialize vectors to store MSEs and stdevs for each machine count
  mse_values_total <- numeric(length(machine_counts))
  sdev_values_total <- numeric(length(machine_counts))
  mse_values_point <- numeric(length(machine_counts))
  sdev_values_point <- numeric(length(machine_counts))
  
  # Loop over each machine count
  for (i in seq_along(machine_counts)) {
    # Print progress
    cat("  Machine count:", machine_counts[i], "\n")

    m <- machine_counts[i]
    ns <- rep(N / m, m)  # Local observations per machine
    
    # Ensure ns are integers
    ns <- floor(ns)
    N_actual <- sum(ns)
    
    # Storage for MSEs
    mse_total_list <- numeric(num_simulations)
    mse_point_list <- numeric(num_simulations)
    
    for (sim in 1:num_simulations) {
      # Print progress
      cat("    Simulation:", sim, "\n")

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

      # plot estimated signal versus true signal
      # plot(x = x_values, y = signal, type = "l", col = "blue", lwd = 2,
      #      main = "Original vs. Reconstructed Signal",
      #      xlab = "x (Time)", ylab = "Amplitude")
      # lines(x = x_values, estimated_signal, col = "red", lwd = 2)
      
      # Calculate MSE over the entire signal
      mse_total <- mean((signal - estimated_signal)^2)
      mse_total_list[sim] <- mse_total

      # Privacy budgets per machine
      pb <- rep(eps, m)
      
      s_point <- s - 1/2
      # Compute parameters
      L_max <- compute_L_max(ns, pb, s_point)
      tau <- sqrt(grid_size) * c_psi + sqrt((2 * s_point + 1) * L_max)
      
      # Estimate signal
      estimated_signal <- federated_estimator(ns, pb, s_point, Y, X, grid_size, max_level, 
        wavelet_family, boundary, S, tau, c_psi)
      
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
    
    # Store results for total MSE
    mse_values_total[i] <- mse_mean_total
    sdev_values_total[i] <- mse_sdev_total
    
    # Store results for point MSE
    mse_values_point[i] <- mse_mean_point
    sdev_values_point[i] <- mse_sdev_point
  }
  
  # Store in results lists
  results_mse_total[[as.character(eps)]] <- list(MSE = mse_values_total, SDEV = sdev_values_total)
  results_mse_point[[as.character(eps)]] <- list(MSE = mse_values_point, SDEV = sdev_values_point)
}

# Function to create LaTeX table and save to file
create_latex_table <- function(results_list, caption_text, label_text, file_name) {
  # Convert results_list to a data frame
  results_df <- data.frame(
    Epsilon = epsilon_values,
    stringsAsFactors = FALSE
  )
  
  # Add columns for each machine count
  for (i in seq_along(machine_counts)) {
    # Print progress
    cat(" Creating table elements for machine count:", machine_counts[i], "\n")

    m <- machine_counts[i]
    
    mse_column <- sapply(results_list, function(x) x$MSE[i])
    sdev_column <- sapply(results_list, function(x) x$SDEV[i])
    
    # Modify the combined strings to include \makecell{} for LaTeX
    combined_strings <- sprintf("\\makecell{%.4f \\\\ (%.4f)}", mse_column, sdev_column)
    
    col_name <- paste0("m=", m)
    results_df[[col_name]] <- combined_strings
  }
  
  # Replace 'Inf' with '\\infty' for LaTeX formatting
  results_df$Epsilon <- as.character(results_df$Epsilon)
  results_df$Epsilon[results_df$Epsilon == "Inf"] <- "$\\infty$"
  
  # Set column names for LaTeX table
  colnames(results_df) <- c("$\\epsilon$", paste0("$m=", machine_counts, "$"))
  
  # Corrected alignment vector
  align_vector <- c("l", rep("c", ncol(results_df)))
  
  # Create LaTeX table
  latex_table <- xtable(results_df, align = align_vector, 
                        caption = caption_text, 
                        label = label_text)
  
  # Capture the LaTeX code
  table_code <- capture.output(print(latex_table, include.rownames = FALSE, sanitize.text.function = identity,
                                     booktabs = TRUE, table.placement = "ht", hline.after = NULL, 
                                     add.to.row = list(pos = list(-1, 0, nrow(results_df)),
                                                       command = c("\\toprule\n", "\\midrule\n", "\\bottomrule\n"))))
  
  # Write the LaTeX code to the specified file
  writeLines(table_code, con = file_name)
  
  # Optionally, print a message indicating that the file has been saved
  cat("LaTeX table saved to", file_name, "\n")
}

# Generate LaTeX table for total MSE and save to file
create_latex_table(results_mse_total, 
                   caption_text = "Mean Squared Error (MSE) over Entire Signal with Standard Deviation Below for Different Privacy Budgets and Machine Counts",
                   label_text = "tab:mse_total",
                   file_name = paste0("N=",N,"mse_total_table.tex"))

# Generate LaTeX table for pointwise MSE and save to file
create_latex_table(results_mse_point, 
                   caption_text = sprintf("Squared Error at Point $t=%.2f$ with Standard Deviation Below for Different Privacy Budgets and Machine Counts", t_point),
                   label_text = "tab:mse_point",
                   file_name = paste0("N=",N,"mse_point_table.tex"))
