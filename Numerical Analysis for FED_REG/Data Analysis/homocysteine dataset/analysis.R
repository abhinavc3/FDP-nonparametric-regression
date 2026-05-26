#' Determine the directory of the currently running script.
#'
#' @return A normalized path to the directory containing the current script, or the current working directory when the script path is unavailable.
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])

  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path[1])))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }

  normalizePath(getwd())
}

script_dir <- get_script_dir()
helpers_dir <- file.path(script_dir, "..", "shared_helpers")
source(file.path(helpers_dir, "wavelet_helper_functions.r"))
source(file.path(helpers_dir, "real_data_plot.r"))
getwd()

# Load necessary libraries
library(foreign)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(haven)

set.seed(2024)


folate_b12_data <- read.xport(file.path(script_dir, "LAB06.XPT"))
thcy_data <- read_xpt(file.path(script_dir, "LAB25.XPT"))

data <- merge(thcy_data, folate_b12_data, by = "SEQN")
# colnames(data)

# Ensure the dataset is in a data frame format
data <- as.data.frame(data)

# Explicitly use dplyr::select to avoid any conflicts
data_filtered <- data %>%
  dplyr::select(SEQN, LBXHCY, LBXFOL, LBXB12) %>% # Select only SEQN, tHcy, Folate, B12
  dplyr::filter(!is.na(LBXHCY), !is.na(LBXFOL), !is.na(LBXB12)) # Remove rows with missing data

# Print the number of rows in the filtered dataset
nrow((data_filtered))

# Rename columns for easier reference
colnames(data_filtered) <- c("ID", "Homocysteine", "Folate", "B12")
head(data_filtered)

# Remove bottom and top 5% quantile from B12 and Folate levels
data_filtered <- data_filtered %>%
  filter(
    B12 > quantile(B12, 0.025), B12 < quantile(B12, 0.975),
    Folate > quantile(Folate, 0.025), Folate < quantile(Folate, 0.975)
  )
summary(data_filtered)
data_filtered_no_outliers <- data_filtered

# Remove rows with missing data
data_filtered_no_outliers <- data_filtered_no_outliers %>%
  filter(!is.na(Homocysteine), !is.na(Folate), !is.na(B12))


# Check the resulting data after outlier removal
summary(data_filtered_no_outliers)
nrow(data_filtered_no_outliers)


# Set parameters for wavelet estimation
grid_size <- 2^15 # Set grid size
S <- 5 # Set smoothness level for the wavelet
wavelet_types <- c("DaubExPhase", "DaubLeAsymm", "Coiflets", "Symmlets", "Haar")
wavelet_family <- wavelet_types[1]
boundary <- "interval"
max_level <- log2(grid_size)


####### PLOT FOR FOLATE VS HOMOCYSTEINE #######

ppX <- data_filtered_no_outliers$Folate
ppY <- data_filtered_no_outliers$Homocysteine
# check for na values
sum(is.na(ppX))
sum(is.na(ppY))

ppY <- log(ppY)
ppX <- log(ppX)

result <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = 3, draw_per_bin = 500)
ppY <- result$ppY
ppX <- result$ppX
#
length(ppX)

pdf(file.path(script_dir, paste0(length(ppX), "folate_vs_homocysteine_plot.pdf")), width = 8, height = 6)
plot_padded_curve(ppX, ppY, s = 2, eps = c(1, 2), S = 5, padding_size_percent = 0.5, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, ylim_range = c(1.2, 2.6), x_lab_name = "log serum Folate (ng/mL)", y_lab_name = "log tHcy (umol/L)")
dev.off()

result <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = 3, draw_per_bin = 250)
ppY <- result$ppY
ppX <- result$ppX
#
length(ppX)

pdf(file.path(script_dir, paste0(length(ppX), "folate_vs_homocysteine_plot.pdf")), width = 8, height = 6)
plot_padded_curve(ppX, ppY, s = 2, eps = c(1, 2), S = 5, padding_size_percent = 0.5, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, ylim_range = c(1.2, 2.6), x_lab_name = "log serum Folate (ng/mL)", y_lab_name = "log tHcy (umol/L)")
dev.off()

result <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = 3, draw_per_bin = 100)
ppY <- result$ppY
ppX <- result$ppX
#
length(ppX)

pdf(file.path(script_dir, paste0(length(ppX), "folate_vs_homocysteine_plot.pdf")), width = 8, height = 6)
plot_padded_curve(ppX, ppY, s = 2, eps = c(1, 2), S = 5, padding_size_percent = 0.5, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, ylim_range = c(1.2, 2.6), x_lab_name = "log serum Folate (ng/mL)", y_lab_name = "log tHcy (umol/L)")
dev.off()

####### PLOT FOR B12 VS HOMOCYSTEINE #######
ppX <- data_filtered_no_outliers$B12
ppY <- data_filtered_no_outliers$Homocysteine
# check for na values
sum(is.na(ppX))
sum(is.na(ppY))

ppY <- log(ppY)
ppX <- log(ppX)

result <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = 3, draw_per_bin = 500)
ppY <- result$ppY
ppX <- result$ppX
#
length(ppX)

pdf(file.path(script_dir, paste0(length(ppX), "b12_vs_homocysteine_plot.pdf")), width = 8, height = 6)
plot_padded_curve(ppX, ppY, s = 2, eps = c(1, 2), S = 5, padding_size_percent = 0.5, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, ylim_range = c(1.2, 2.6), x_lab_name = "log serum B12 (pg/mL)", y_lab_name = "log tHcy (umol/L)")
dev.off()

result <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = 3, draw_per_bin = 250)
ppY <- result$ppY
ppX <- result$ppX
#
length(ppX)

pdf(file.path(script_dir, paste0(length(ppX), "b12_vs_homocysteine_plot.pdf")), width = 8, height = 6)
plot_padded_curve(ppX, ppY, s = 2, eps = c(1, 2), S = 5, padding_size_percent = 0.5, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, ylim_range = c(1.2, 2.6), x_lab_name = "log serum B12 (pg/mL)", y_lab_name = "log tHcy (umol/L)")
dev.off()

result <- subsample_uniform_ppX_ppY(ppX, ppY, num_bins = 3, draw_per_bin = 100)
ppY <- result$ppY
ppX <- result$ppX
#
length(ppX)

pdf(file.path(script_dir, paste0(length(ppX), "b12_vs_homocysteine_plot.pdf")), width = 8, height = 6)
plot_padded_curve(ppX, ppY, s = 2, eps = c(1, 2), S = 5, padding_size_percent = 0.5, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, ylim_range = c(1.2, 2.6), x_lab_name = "log serum B12 (pg/mL)", y_lab_name = "log tHcy (umol/L)")
dev.off()
