#setwd("/Users/lassev/Dropbox (Penn)/Distributed comm. and privacy constraints/Federated Learning for Nonparametric Regression/Simulations/expiratory/")
source("../wavelet_helper_functions.r")
source("../real_data_plot.r")
getwd()

# Load necessary libraries
library(foreign)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(haven)

# Set seed
set.seed(2024)

# Load the data
data_fev <- read_xpt("SPX_F.XPT")
data_demographics <- read_xpt("DEMO_F.XPT")
# Merge the data using the common identifier 'SEQN'
merged_data <- inner_join(data_fev, data_demographics, by = "SEQN")
#colnames(merged_data)

# Extract variables for plotting (adjust variable names as needed)
age <- merged_data$RIDAGEYR
fev1 <- merged_data$SPXNFEV1/1000  # Convert to liters
gender <- merged_data$RIAGENDR

# Create a data frame for plotting
plot_data <- data.frame(age, fev1, gender)

# Remove missing values
plot_data <- na.omit(plot_data)

# Filter out unrealistic values for age and FEV1
plot_data <- plot_data %>% filter(age >= 6, age <= 80, fev1 > 0)

# Stratified Wavelet Estimation by Gender
plot_data_male <- plot_data %>% filter(gender == 1)
plot_data_female <- plot_data %>% filter(gender == 2)

# Total number of observations
nrow(plot_data)
nrow(plot_data_male)
nrow(plot_data_female)


# Set parameters for wavelet estimation
grid_size <- 2^15  # Set grid size
S <- 8  # Set smoothness level for the wavelet
wavelet_types <- c("DaubExPhase", "DaubLeAsymm", "Coiflets", "Symmlets", "Haar")
wavelet_family <- wavelet_types[1]  
boundary <- "interval"
max_level <- log2(grid_size)

# Wavelet Estimation Preparation
ppX <- plot_data$age
ppY <- plot_data$fev1

#plot_padded_curve(ppX, ppY, s=2, eps=c(0.5,1,2), S=5, padding_size_percent = 0.05, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = FALSE, ylim_range = c(0.5, 5), symbol_spacing = 2^11, x_lab_name = "Age (in years)", y_lab_name = "base-10 log of FEV1 (in liters)")

# Wavelet Estimation Preparation
ppX1 <- plot_data_female$age
ppY1 <- plot_data_female$fev1
ppX2 <- plot_data_male$age
ppY2 <- plot_data_male$fev1
length(ppX1)
length(ppX2)

pdf("female_male_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "female", dataset2_name = "male")
dev.off()

# INCLUDING SMOKING STATUS

# Load the data
data_fev <- read_xpt("SPX_F.XPT")
data_demographics <- read_xpt("DEMO_F.XPT")
data_smoking <- read_xpt("SMQ_F.XPT")
# Merge the data using the common identifier 'SEQN'
merged_data <- inner_join(data_fev, data_demographics, by = "SEQN")
merged_data <- inner_join(merged_data, data_smoking, by = "SEQN") 

age <- merged_data$RIDAGEYR
fev1 <- merged_data$SPXNFEV1/1000  # Convert to liters
gender <- merged_data$RIAGENDR 
#colnames((data_smoking))
smoker_at_least_100_lifetime <- merged_data$SMQ020
smoker_now_smoke <- merged_data$SMQ040

plot_data_smoke_lt <- data.frame(age, fev1, gender, smoker_at_least_100_lifetime)

# Remove missing values
plot_data_smoke_lt <- na.omit(plot_data_smoke_lt)

plot_data_smoke_lt_smoker <- plot_data_smoke_lt %>% filter(smoker_at_least_100_lifetime == 1)
plot_data_smoke_lt_nonsmoker <- plot_data_smoke_lt %>% filter(smoker_at_least_100_lifetime == 2)

# Seperate the responses by gender
plot_data_smoke_lt_smoker_men <- plot_data_smoke_lt_smoker %>% filter(gender == 1)
plot_data_smoke_lt_smoker_women <- plot_data_smoke_lt_smoker %>% filter(gender == 2)
plot_data_smoke_lt_nonsmoker_men <- plot_data_smoke_lt_nonsmoker %>% filter(gender == 1)
plot_data_smoke_lt_nonsmoker_women <- plot_data_smoke_lt_nonsmoker %>% filter(gender == 2)

ppX1 <- plot_data_smoke_lt_smoker_women$age
ppY1 <- plot_data_smoke_lt_smoker_women$fev1
length(ppX1)
#plot_padded_curve(ppX1, ppY1, s=2, eps=c(5,10), S=5, padding_size_percent = 0.05, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11)

ppX2 <- plot_data_smoke_lt_nonsmoker_women$age
ppY2 <- plot_data_smoke_lt_nonsmoker_women$fev1
length(ppX2)
#plot_padded_curve(ppX2, ppY2, s=2, eps=c(5,10), S=5, padding_size_percent = 0.05, grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11)

pdf("smoking_female_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "smoking", dataset2_name = "non-smoking")
dev.off()


ppX1 <- plot_data_smoke_lt_smoker_men$age
ppY1 <- plot_data_smoke_lt_smoker_men$fev1

ppX2 <- plot_data_smoke_lt_nonsmoker_men$age
ppY2 <- plot_data_smoke_lt_nonsmoker_men$fev1

pdf("smoking_male_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "smoking", dataset2_name = "non-smoking")
dev.off()


# Create a data frame for plotting
plot_data_smoke_lt <- data.frame(age, fev1, gender, smoker_at_least_100_lifetime)
length(age)
length(smoker_at_least_100_lifetime)
# Remove missing values
plot_data_smoke_lt <- na.omit(plot_data_smoke_lt)
dim(plot_data_smoke_lt)
# Separate the responses by smoking status
plot_data_smoke_lt_smoker <- plot_data_smoke_lt %>% filter(smoker_at_least_100_lifetime == 1)
plot_data_smoke_lt_nonsmoker <- plot_data_smoke_lt %>% filter(smoker_at_least_100_lifetime == 2)

plot_data_smoke_now <- data.frame(age, fev1, gender, smoker_now_smoke)
plot_data_smoke_now <- na.omit(plot_data_smoke_now)
plot_data_smoke_now_smoker <- plot_data_smoke_now %>% filter(smoker_now_smoke == 1)
plot_data_smoke_now_nonsmoker <- plot_data_smoke_now %>% filter(smoker_now_smoke == 2)

ppX1 <- plot_data_smoke_now_smoker$age
ppY1 <- plot_data_smoke_now_smoker$fev1

ppX2 <- plot_data_smoke_now_nonsmoker$age
ppY2 <- plot_data_smoke_now_nonsmoker$fev1

length(ppX1)
length(ppX2)

pdf("smoking_nonsmoking_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "every day", dataset2_name = "some days", ylim_range = c(2,5))
dev.off()


### SMOKING VERSUS NON-SMOKING

plot_data_smoke_now <- data.frame(age, fev1, gender, smoker_now_smoke)
plot_data_smoke_now <- na.omit(plot_data_smoke_now)
plot_data_smoke_now_smoker <- plot_data_smoke_now %>% filter(smoker_now_smoke == 1)
plot_data_smoke_now_nonsmoker <- plot_data_smoke_now %>% filter(smoker_now_smoke == 3)

ppX1 <- plot_data_smoke_now_smoker$age
ppY1 <- plot_data_smoke_now_smoker$fev1

ppX2 <- plot_data_smoke_now_nonsmoker$age
ppY2 <- plot_data_smoke_now_nonsmoker$fev1

length(ppX1)
length(ppX2)

pdf("smoking_fully_nonsmoking_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "every day", dataset2_name = "never", ylim_range = c(2,5))
dev.off()

# SEGMENT BY GENDER

plot_data_smoke_now_smoker_male <- plot_data_smoke_now_smoker %>% filter(gender == 1)
plot_data_smoke_now_nonsmoker_male <- plot_data_smoke_now_nonsmoker %>% filter(gender == 1)

ppX1 <- plot_data_smoke_now_smoker_male$age
ppY1 <- plot_data_smoke_now_smoker_male$fev1

ppX2 <- plot_data_smoke_now_nonsmoker_male$age
ppY2 <- plot_data_smoke_now_nonsmoker_male$fev1

length(ppX1)
length(ppX2)

pdf("smoking_male_fully_nonsmoking_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "every day", dataset2_name = "never", ylim_range = c(2,5))
dev.off()


# SEGMENT BY GENDER - FEMALE

plot_data_smoke_now_smoker_female <- plot_data_smoke_now_smoker %>% filter(gender == 2)
plot_data_smoke_now_nonsmoker_female <- plot_data_smoke_now_nonsmoker %>% filter(gender == 2)

ppX1 <- plot_data_smoke_now_smoker_male$age
ppY1 <- plot_data_smoke_now_smoker_male$fev1

ppX2 <- plot_data_smoke_now_nonsmoker_male$age
ppY2 <- plot_data_smoke_now_nonsmoker_male$fev1

length(ppX1)
length(ppX2)

pdf("smoking_female_fully_nonsmoking_fev1_plot.pdf", width = 8, height = 6)
plot_combined_curves_two_datasets(ppX1, ppY1, ppX2, ppY2, s = 2, epsilons = c(1), S = 5, tau = 6, padding_size_percent = 0.05, x_lab_name = "Age (years)", y_lab_name = "FEV1 (liters)", grid_size = 2^15, wavelet_family = "DaubExPhase", boundary = "interval", include_scatter = TRUE, symbol_spacing = 2^11, dataset1_name = "every day", dataset2_name = "never", ylim_range = c(2,5))
dev.off()

