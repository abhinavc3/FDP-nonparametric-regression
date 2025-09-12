library(dplyr)
library(ggplot2)
library(tidyr)
library(Metrics)
#require(mvmesh)

source('method-functions-fastcov.R')

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

# Disease status count by source
ggplot(combined_data, aes(x = factor(disease), fill = source)) +
  geom_bar(position = "dodge") +
  labs(title = "Disease Status Distribution by Source", x = "Disease Status", y = "Count") +
  theme_minimal()


# Age Distribution by Source
ggplot(combined_data, aes(x = age, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Age Distribution by Source", x = "Age", y = "Frequency") +
  theme_minimal()

# Cholesterol Levels by Source
ggplot(combined_data, aes(x = chol, fill = source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Cholesterol Levels by Source", x = "Cholesterol", y = "Density") +
  theme_minimal()

# Resting Blood Pressure by Source
ggplot(combined_data, aes(x = trestbps, fill = source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Resting Blood Pressure by Source", x = "Resting Blood Pressure", y = "Density") +
  theme_minimal()

# Max Heart Rate by Disease Status and Source
ggplot(combined_data, aes(x = factor(disease), y = thalach, fill = source)) +
  geom_boxplot() +
  labs(title = "Max Heart Rate by Disease Status and Source", x = "Disease Status", y = "Max Heart Rate") +
  theme_minimal()

f1score <- function(pred, tru){
  
  if(mean(tru) > 0.5){
    pred <- 1 - pred
    tru <- 1 - tru
  }
  
  tp <- sum(pred*tru)
  fp <- sum(pred*(1-tru))
  fn <- sum((1-pred)*tru)
  return(2*tp/(2*tp+fp+fn))
}

scaler <- function(xmat){
  d <- ncol(xmat)
  for(j in 1:d){
    
    jrange <- max(new_data1[,j]) - min(new_data1[,j])
    if(jrange == 0){
      xmat[,j] = jrange * xmat[,j]
    }else{
      xmat[,j] <- 0.5*(xmat[,j]-min(new_data1[,j]))/jrange}
  }
  return(xmat)
}

twoclass <- function(y){
  y = c(y)
  y[y>1] = 1
  return(y)
}

varset <- seq(1:15)[-c(5,6,11,12,13)]

new_data <- (combined_data[,varset])

new_data0 <- data.frame(new_data$age,
                        new_data$sex,
                        #new_data$cp,#+0.2*runif(length(new_data$cp)),
                        rowMeans(diag(4)[new_data$cp,1:2]), 
                        #new_data$cp,
                        new_data$exang,
                        new_data$thalach,
                        new_data$oldpeak,
                        new_data$trestbps,
                        #new_data$restecg,
                        #(diag(3)[new_data$restecg+1,2]),
                        new_data$disease,
                        new_data$source)

colnames(new_data0) <- c("age",
                         "sex",
                         "cp",
                         #"cp1", "cp2", "cp3", "cp4",
                         "exang",
                         "thalach",
                         "oldpeak",
                         "trestbps",
                         #"restecg",
                         #"recg0", "recg1", "recg2", 
                         "disease", "source")

npreds <- ncol(new_data0)-2

new_data1 <- na.omit(new_data0)

servers <- c("Cleveland", "Hungarian", "Switzerland", "VA")

target_idx <- 2

target_data <- new_data1%>%filter(new_data1$source == servers[target_idx])

sources <- servers[-target_idx]

logist_acc <- function(){
  testrows <- sample(seq(1:nrow(target_data)), 100)
  
  test_data <- target_data[testrows,]
  
  train_data <- target_data[-testrows,]
  
  for(i in 1:length(sources)){
    train_data <- bind_rows(train_data,
                            new_data1%>%filter(new_data1$source == sources[i]))
  }
  
  #testrows <- sample(seq(1:303),100)
  logist <- glm(disease~., data = train_data, family = binomial)
  
  logis_pred <- predict(logist, newdata = test_data)
  ac_res <- function(th){
    mean(ifelse(logis_pred > th,1,0) == (new_data1[testrows,])$disease)}
  thseq <- seq(-1,1,0.1)
  ypred <- thseq
  for(i in 1:length(thseq)){
    ypred[i] <- ac_res(thseq[i])
  }
  return(max(ypred))
}

nservs <- length(servers)
#wmat <- UnitSimplex(nservs, 12)$V

wmat <- read.csv('wmat.csv', header = T)

wmat <- as.matrix(wmat[,-1])
dimnames(wmat) = NULL

heart_accuracy <- function(eps1){
  
  testrows <- sample(seq(1:nrow(target_data)), max(150, nrow(target_data)/2))
  
  train_targetX <- data.matrix(target_data[-testrows, seq(1:npreds)])
  train_targetY <- data.matrix(target_data[-testrows,]$disease)
  
  dimnames(train_targetX) <- NULL
  dimnames(train_targetY) <- NULL
  
  train_target <- list( X = scaler(train_targetX), Y = twoclass(train_targetY))
  
  source1_data <- new_data1%>%filter(new_data1$source == sources[1])
  
  train_source1X <- data.matrix(source1_data[, seq(1:npreds)])
  train_source1Y <- data.matrix(source1_data$disease)
  dimnames(train_source1X) <- NULL
  dimnames(train_source1Y) <- NULL
  
  source2_data <- new_data1%>%filter(new_data1$source == sources[2])
  train_source2X <- data.matrix(source2_data[, seq(1:npreds)])
  
  train_source2Y <- data.matrix(source2_data$disease)
  dimnames(train_source2X) <- NULL
  dimnames(train_source2Y) <- NULL
  
  source3_data <- new_data1%>%filter(new_data1$source == sources[3])
  train_source3X <- data.matrix(source3_data[, seq(1:npreds)])
  
  train_source3Y <- data.matrix(source3_data$disease)
  dimnames(train_source3X) <- NULL
  dimnames(train_source3Y) <- NULL
  
  train_source <- list( list( X = scaler(train_source1X), Y = twoclass(train_source1Y)),
                        list( X = scaler(train_source2X), Y = twoclass(train_source2Y)),
                        list( X = scaler(train_source3X), Y = twoclass(train_source3Y)))
  
  # Generate test data
  
  testX <- data.matrix(target_data[testrows, seq(1:npreds)])
  testY <- data.matrix(target_data[testrows,]$disease)
  
  dimnames(testX) <- NULL
  dimnames(testY) <- NULL
  
  test <- list( X = scaler(testX), Y = twoclass(testY))
  
  samp_sizes <- c(nrow(train_targetX), nrow(train_source1X), nrow(train_source2X), nrow(train_source3X))
  
  samp_props <- samp_sizes/sum(samp_sizes)
  
  kernel_f = triangular_kernel
  
  pred_ada <- p_kern_adapt(test , train_target = train_target, train_source = train_source, kernel_f = kernel_f, eps = eps1, wmat = wmat)
  pred_ada_tar <- p_kern_adapt(test , train_target = train_target, train_source = train_source, kernel_f = kernel_f, target_only = T, eps = eps1)
  #pred_ada_eq <- p_kern_adapt(test , train_target = train_target, train_source = train_source, kernel_f = kernel_f, wvec = rep(1/4, 4), eps = eps1)
  pred_ada_samp <- p_kern_adapt(test , train_target = train_target, train_source = train_source, kernel_f = kernel_f, wvec = samp_props, eps = eps1)
  pred_ada_homog <- p_kern_adapt(test , train_target = train_target, train_source = train_source, kernel_f = kernel_f, hmognus = T, eps = eps1)
  
  
  mean_res <- rbind(mean(pred_ada[[2]] == (test$Y)),
                mean(pred_ada_tar[[2]] == (test$Y)),
                #mean(pred_ada_eq[[2]] == (test$Y)),
                mean(pred_ada_samp[[2]] == (test$Y)),
                mean(pred_ada_homog[[2]] == (test$Y))
                )
  
  f1_res <- rbind(f1score(pred_ada[[2]],test$Y),
                    f1score(pred_ada_tar[[2]],test$Y),
                    #f1score(pred_ada_eq[[2]],test$Y),
                    f1score(pred_ada_samp[[2]],test$Y),
                    f1score(pred_ada_homog[[2]],test$Y)
                    )
  
  
  auc_res <- rbind(auc(test$Y,pred_ada[[1]]),
                   auc(test$Y,pred_ada_tar[[1]]),
                   #auc(test$Y,pred_ada_eq[[1]]),
                   auc(test$Y,pred_ada_samp[[1]]),
                   auc(test$Y,pred_ada_homog[[1]])
  )
  return(data.frame(accuracy = mean_res, f1scores = f1_res, auc = auc_res))
  
}
function_summary <- function(x) {
  mean = mean(x)
  sd = sd(x) 
  return(c(
    mean = mean,
    sd = sd,
    lower_quantile = mean - 1.96 * sd / sqrt(length(x)),
    upper_quantile = mean + 1.96 * sd / sqrt(length(x))
  ))
}
sim_summary <- function(sim_result, t = 3) {
  cbind(sim_result[, 1:t], t(apply(sim_result[,-(1:t)], 1, function_summary)))
}

eps_grid = c(1,2,4,8,16)
num_rep = 200
#methods <- c("AdaptAll", "AdaptTar", "AdaptEq", "AdaptSamp", "AdaptHomog")
methods <- c("AdaptAll", "AdaptTar", "AdaptSamp", "AdaptHomog")
acc_results <- data.frame(matrix(nrow = length(eps_grid) * length(methods), ncol = num_rep + 2))
f1_results <- data.frame(matrix(nrow = length(eps_grid) * length(methods), ncol = num_rep + 2))
auc_results <- data.frame(matrix(nrow = length(eps_grid) * length(methods), ncol = num_rep + 2))

index = 1 
for(eps in eps_grid){
  output = replicate(num_rep, heart_accuracy(eps))
  for (i in 1:length(methods)) {
    acc_results[index,1] <- eps
    acc_results[index,2] <- methods[i]
    acc_results[index, -(1:2)] <- sapply(output['accuracy',], "[", i)
    
    auc_results[index,1] <- eps
    auc_results[index,2] <- methods[i]
    auc_results[index, -(1:2)] <- sapply(output['auc',], "[", i)
    
    f1_results[index,1] <- eps
    f1_results[index,2] <- methods[i]
    f1_results[index, -(1:2)] <- sapply(output['f1scores',], "[", i)
    
    
    index <- index + 1
  }
  print(eps)
}
final_acc_results =sim_summary(acc_results, t = 2)
colnames(final_acc_results)[1:2] <- c("eps", "method")

final_auc_results =sim_summary(auc_results, t = 2)
colnames(final_auc_results)[1:2] <- c("eps", "method")

final_f1scores_results =sim_summary(f1_results, t = 2)
colnames(final_f1scores_results)[1:2] <- c("eps", "method")

# Plot with whiskers
# Plot with colored whiskers
ggplot(final_acc_results, aes(x = eps, y = mean, color = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile, color = method), width = 0.1) +  # Mapping color aesthetic
  labs(x = expression(epsilon), y = "Mean Accuracy", color = "Method") +
  theme_minimal()

## AUc plot
ggplot(final_auc_results, aes(x = eps, y = mean, color = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile, color = method), width = 0.1) +  # Mapping color aesthetic
  labs(x = expression(epsilon), y = "Mean AUC", color = "Method") +
  theme_minimal()

## f1scores plot
ggplot(final_f1scores_results, aes(x = eps, y = mean, color = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile, color = method), width = 0.1) +  # Mapping color aesthetic
  labs(x = expression(epsilon), y = "Mean f1 Scores", color = "Method") +
  theme_minimal()

library(ggplot2)
library(patchwork)

color_palette <- c(
  "AdaptTar" = "#E69F00",  
  #"AdaptEq" = "#56B4E9",
  "AdaptAll" = "#9900CC",  
  "AdaptSamp" = "#008080",  
  "AdaptHomog" = "#56B4E9"  
)

# Custom theme for bold and larger fonts
custom_theme <- function() {
  theme_minimal() +
    theme(text = element_text(size = 10),  # Default text size and style
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "none")
}

# Plot for accuracy
p1 <- ggplot(final_acc_results, aes(x = eps, y = mean, color = method, shape = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile), width = 0.01) +
  labs(x = expression(epsilon), y = "Mean Accuracy") +
  scale_color_manual(values = color_palette, name = "Method") +  
  scale_shape_manual(values = c(18, 17, 21, 19, 22), name = "Method") +
  custom_theme()

# Plot for AUC
p2 <- ggplot(final_auc_results, aes(x = eps, y = mean, color = method, shape = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile), width = 0.01) +
  labs(x = expression(epsilon), y = "Mean AUC") +
  scale_color_manual(values = color_palette, name = "Method") +  
  scale_shape_manual(values = c(18, 17, 21, 19, 20), name = "Method") +
  custom_theme()

# Plot for f1scores
p3 <- ggplot(final_f1scores_results, aes(x = eps, y = mean, color = method, shape = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower_quantile, ymax = upper_quantile), width = 0.01) +
  labs(x = expression(epsilon), y = "Mean F1 Score") +
  scale_color_manual(values = color_palette, name = "Method") +  
  scale_shape_manual(values = c(18, 17, 21, 19, 22), name = "Method") +
  custom_theme()


# Combine the plots side by side and add a single legend
combined_plot <- p1 + p3 + plot_layout(guides = 'collect') +
  theme(legend.position = "right")

# Print the combined plot
print(combined_plot)


#mean_logist <- replicate(100, logist_acc())

#mean(mean_logist)
#res_df





