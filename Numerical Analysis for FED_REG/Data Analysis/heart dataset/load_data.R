# --- Libraries and Data Loading ---
library(dplyr)
library(ggplot2)
library(tidyr)
library(Metrics)
library(GGally)
library(RColorBrewer)
library(wavethresh)  # for wavelet routines

# Load the heart disease datasets
clev_data <- read.csv('heart+disease/processed.cleveland.data', header = FALSE)
swiss_data <- read.csv('heart+disease/processed.switzerland.data', header = FALSE)
hung_data <- read.csv('heart+disease/processed.hungarian.data', header = FALSE)
va_data <- read.csv('heart+disease/processed.va.data', header = FALSE)

# Add a source column and combine datasets
clev_data <- clev_data %>% mutate(source = "Cleveland")
swiss_data <- swiss_data %>% mutate(source = "Switzerland")
hung_data <- hung_data %>% mutate(source = "Hungarian")
va_data <- va_data %>% mutate(source = "VA")
combined_data <- rbind(clev_data, swiss_data, hung_data, va_data)

# Assign column names
column_names <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
                  "thalach", "exang", "oldpeak", "slope", "ca", "thal", "disease", "source")
colnames(combined_data) <- column_names

# Recode disease: presence if disease > 1
combined_data$disease <- ifelse(combined_data$disease > 1, 1, combined_data$disease)

# Replace "?" with NA and convert columns to numeric (except source)
combined_data <- combined_data %>% 
  mutate(across(where(is.character), ~na_if(., "?"))) %>%
  mutate(across(-source, ~as.numeric(as.character(.))))

# --- Subsample to Obtain Uniform X Distribution ---
# In this example, we reconstruct the relationship between age and thalach.
# (You can filter out any missing values first.)
data_filtered <- combined_data %>% 
  filter(!is.na(age) & !is.na(thalach))

