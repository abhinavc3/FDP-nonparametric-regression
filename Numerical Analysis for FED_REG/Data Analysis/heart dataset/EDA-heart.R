library(dplyr)
library(ggplot2)
library(tidyr)
library(Metrics)
library(GGally)
#require(mvmesh)


clev_data <- read.csv('heart+disease/processed.cleveland.data', header = F)
head(clev_data)

swiss_data <- read.csv('heart+disease/processed.switzerland.data', header = F)
head(swiss_data)

hung_data <- read.csv('heart+disease/processed.hungarian.data', header = F)
head(hung_data)

va_data <- read.csv('heart+disease/processed.va.data', header = F)
head(va_data)



# Adding a source column to each data frame
clev_data <- clev_data %>% mutate(source = "Cleveland")
swiss_data <- swiss_data %>% mutate(source = "Switzerland")
hung_data <- hung_data %>% mutate(source = "Hungarian")
va_data <- va_data %>% mutate(source = "VA")

# Combine the data frames
combined_data <- rbind(clev_data, swiss_data, hung_data, va_data)

# View the first few rows of the combined data frame
head(combined_data)

# Column names for the data
column_names <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
                  "thalach", "exang", "oldpeak", "slope", "ca", "thal", "disease", "source")

# Assign these names to the combined_data frame
colnames(combined_data) <- column_names
# Assuming combined_data is your data frame and 'disease' is the column
combined_data$disease <- ifelse(combined_data$disease > 1, 1, combined_data$disease)

head(combined_data)



combined_data <- combined_data %>%
  mutate(across(where(is.character), ~na_if(., "?")))

# For the 'chol' column and any other specific cases where "0" is considered a missing value:
combined_data$chol <- as.numeric(na_if(combined_data$chol, "0"))

str(combined_data)


# Assuming combined_data is your DataFrame
# Convert all columns to numeric except 'source'
combined_data <- combined_data %>%
  mutate(across(-source, ~as.numeric(as.character(.))))

# Check the structure
str(combined_data)

# View summary statistics
summary(combined_data)

# Visualize the missing data proportion
missing_data <- combined_data %>%
  summarise(across(everything(), ~mean(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "na_proportion")

ggplot(missing_data, aes(x = variable, y = na_proportion, fill = variable)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Variable", y = "Proportion Missing", title = "Missing Data Proportions by Variable")

# Calculate missing data proportions grouped by 'source'
missing_data <- combined_data %>%
  group_by(source) %>%  # Group data by 'source'
  summarise(across(everything(), ~mean(is.na(.)), .names = "na_{.col}")) %>%
  pivot_longer(cols = starts_with("na_"), names_to = "variable", values_to = "na_proportion") %>%
  mutate(variable = sub("na_", "", variable))  # Clean column names for plotting

# Plot missing data proportions faceted by 'source'
ggplot(missing_data, aes(x = variable, y = na_proportion, fill = variable)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ source) +  # Create a separate plot for each source
  theme_minimal() +
  labs(x = "Variable", y = "Proportion Missing", title = "Missing Data Proportions by Variable and Source") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Plot Age vs Cholesterol
ggplot(combined_data, aes(x = age, y = chol, color = source)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Age vs Cholesterol", x = "Age", y = "Cholesterol") +
  theme(legend.position = "bottom")

# Plot Age vs Blood Pressure
ggplot(combined_data, aes(x = age, y = trestbps, color = source)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Age vs Blood Pressure", x = "Age", y = "Resting Blood Pressure") +
  theme(legend.position = "bottom")

# Plot Age vs Thalach
ggplot(combined_data, aes(x = age, y = thalach, color = source)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Age vs Thalach", x = "Age", y = "Max Heart Rate") +
  theme(legend.position = "bottom")

continuous_vars <- combined_data[c('age', 'trestbps', 'chol', 'thalach', 'oldpeak')]

# Pair plot of continuous variables
ggpairs(continuous_vars, 
        lower = list(continuous = wrap("points", alpha = 0.5)),
        title = "Pairwise Scatterplots of Continuous Features")

# Age Distribution by Source
ggplot(combined_data, aes(x = age, fill = source)) +
  geom_density(alpha = 0.5) +
  labs( x = "Age", y = "Frequency") +
  theme_minimal()


# Ensure the data is numeric
combined_data <- combined_data %>%
  mutate(age = as.numeric(age),
         thalach = as.numeric(thalach))

# 1. Plot thalach vs age with smooth curves for each source (faceted)
p1 <- ggplot(combined_data, aes(x = age, y = thalach, color = source)) +
  geom_point(alpha = 0.5) +  # Scatter plot with transparency
  geom_smooth(method = "loess", se = FALSE) +  # Smoothed LOESS curve
  facet_wrap(~source) +  # Creates separate plots for each site
  ylim(0, 250) +  # Set Y-axis range
  theme_minimal() +
  labs(title = "Thalach vs Age Across Different Sites",
       x = "Age",
       y = "Max Heart Rate Achieved (Thalach)") +
  theme(legend.position = "none")  # Remove legend since it's faceted

# 2. Combined plot with all data and one smoothed curve
p2 <- ggplot(combined_data, aes(x = age, y = thalach)) +
  geom_point(alpha = 0.3, color = "black") +  # Scatter plot with lower transparency
  geom_smooth(method = "loess", color = "blue", se = FALSE) +  # One smooth LOESS curve
  ylim(0, 250) +  # Set Y-axis range
  theme_minimal() +
  labs(title = "Thalach vs Age (All Sites Combined)",
       x = "Age",
       y = "Max Heart Rate Achieved (Thalach)")

# Print both plots
print(p1)  # Faceted plot
print(p2)  # Combined plot
