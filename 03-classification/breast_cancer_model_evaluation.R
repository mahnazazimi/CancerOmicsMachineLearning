# Breast Cancer Model Evaluation

library(caret)
library(ggplot2)
library(e1071)
library(pROC)
library(rpart)
library(rpart.plot)

data <- read.csv("breast-cancer.csv")

# Use sample IDs as row names and remove the ID column
rownames(data) <- data[, 1]
data <- data[, -1]

dim(data)
table(data$diagnosis)

# Convert diagnosis to a categorical outcome
# B = Benign, M = Malignant
data$diagnosis <- factor(
  data$diagnosis,
  levels = c("B", "M"),
  labels = c("Benign", "Malignant")
)



# Stratified train/test split
train_index <- createDataPartition(
  data$diagnosis,
  p = 0.8,
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

table(train_data$diagnosis)
table(test_data$diagnosis)


# Logistic Regression
lr_model <- glm(
  diagnosis ~ .,
  data = train_data,
  family = binomial
)

# Probability of malignant diagnosis
lr_probability <- predict(
  lr_model,
  newdata = test_data,
  type = "response"
)

hist(
  lr_probability,
  main = "Logistic Regression Probabilities",
  xlab = "Predicted probability of malignancy"
)

# Convert probabilities to class predictions
lr_pred <- ifelse(
  lr_probability > 0.5,
  "Malignant",
  "Benign"
)

lr_pred <- factor(
  lr_pred,
  levels = levels(data$diagnosis)
)


# Logistic Regression Evaluation
lr_conf_matrix <- confusionMatrix(
  lr_pred,
  test_data$diagnosis,
  positive = "Malignant"
)

lr_conf_matrix

conf_matrix_data <- as.data.frame(
  as.table(lr_conf_matrix$table)
)

ggplot(
  conf_matrix_data,
  aes(Prediction, Reference, fill = Freq)
) +
  geom_tile() +
  geom_text(
    aes(label = Freq),
    color = "black",
    size = 6
  ) +
  scale_fill_gradient(
    low = "khaki",
    high = "maroon"
  ) +
  labs(
    x = "Predicted",
    y = "Actual",
    title = "Logistic Regression Confusion Matrix"
  ) +
  theme_minimal()


# Logistic Regression Feature Importance
coefficients <- summary(lr_model)$coefficients

feature_df <- data.frame(
  Feature = rownames(coefficients),
  Coefficient = coefficients[, "Estimate"]
)

# Intercept is not a predictor feature
feature_df <- feature_df[
  feature_df$Feature != "(Intercept)",
]

feature_df$Importance <- abs(
  feature_df$Coefficient
)

feature_df <- feature_df[
  order(feature_df$Importance, decreasing = TRUE),
]

rownames(feature_df) <- NULL

View(feature_df)

ggplot(
  feature_df,
  aes(
    x = reorder(Feature, Importance),
    y = Importance,
    fill = Importance
  )
) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "gold3",
    high = "maroon4"
  ) +
  coord_flip() +
  labs(
    title = "Feature Importance in Logistic Regression",
    x = "Features",
    y = "Absolute Coefficient"
  )


# Support Vector Machine
svm_model <- svm(
  diagnosis ~ .,
  data = train_data,
  probability = TRUE
)

svm_pred <- predict(
  svm_model,
  newdata = test_data,
  probability = TRUE
)

# Extract class probabilities
svm_probabilities <- attr(
  svm_pred,
  "probabilities"
)

svm_probabilities <- as.data.frame(
  svm_probabilities
)

# Probability of malignant diagnosis
svm_malignant_probability <- svm_probabilities[
  , "Malignant"
]


# SVM Evaluation: ROC and AUC
svm_roc <- roc(
  response = test_data$diagnosis,
  predictor = svm_malignant_probability,
  levels = c("Benign", "Malignant")
)

svm_auc <- auc(svm_roc)
svm_auc

plot(
  svm_roc,
  col = "purple4",
  lwd = 3,
  main = "SVM ROC Curve"
)

text(
  x = 0.8,
  y = 0.8,
  label = paste0(
    "AUC: ",
    round(svm_auc, 3)
  )
)


# Sensitivity, Specificity and Best Threshold
sensitivity <- svm_roc$sensitivities
specificity <- svm_roc$specificities
thresholds  <- svm_roc$thresholds

roc_values <- data.frame(
  Threshold = thresholds,
  Sensitivity = sensitivity,
  Specificity = specificity
)

View(roc_values)

# Find threshold with the highest combined
# sensitivity and specificity
best_index <- which.max(
  sensitivity + specificity
)

best_sensitivity <- sensitivity[best_index]
best_specificity <- specificity[best_index]
best_threshold <- thresholds[best_index]

best_threshold

points(
  best_specificity,
  best_sensitivity,
  pch = 19,
  col = "tomato"
)

hist(
  svm_malignant_probability,
  main = "SVM Predicted Probabilities",
  xlab = "Probability of malignancy"
)

abline(
  v = best_threshold,
  col = "darkred"
)

svm_pred_best <- ifelse(
  svm_malignant_probability >= best_threshold,
  "Malignant",
  "Benign"
)

svm_pred_best <- factor(
  svm_pred_best,
  levels = levels(data$diagnosis)
)

svm_conf_matrix <- confusionMatrix(
  svm_pred_best,
  test_data$diagnosis,
  positive = "Malignant"
)

svm_conf_matrix


# SVM Feature Importance
svm_caret_model <- train(
  diagnosis ~ .,
  data = train_data,
  method = "svmRadial"
)

svm_importance <- varImp(
  svm_caret_model,
  scale = TRUE
)

svm_feature_df <- data.frame(
  Feature = rownames(svm_importance$importance),
  Importance = svm_importance$importance[, 1]
)

rownames(svm_feature_df) <- NULL

View(svm_feature_df)

ggplot(
  svm_feature_df,
  aes(
    x = reorder(Feature, Importance),
    y = Importance,
    fill = Importance
  )
) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "gold3",
    high = "maroon4"
  ) +
  coord_flip() +
  labs(
    title = "Feature Importance in SVM",
    x = "Features",
    y = "Importance"
  )


# Decision Tree
dt_model <- rpart(
  diagnosis ~ .,
  data = train_data,
  method = "class"
)

# Class prediction
dt_pred <- predict(
  dt_model,
  newdata = test_data,
  type = "class"
)

# Probability prediction
dt_probability <- predict(
  dt_model,
  newdata = test_data,
  type = "prob"
)


# Decision Tree Evaluation
dt_conf_matrix <- confusionMatrix(
  dt_pred,
  test_data$diagnosis,
  positive = "Malignant"
)

dt_conf_matrix


# Decision Tree Visualization
rpart.plot(dt_model)


# Decision Tree Feature Importance
dt_caret_model <- train(
  diagnosis ~ .,
  data = train_data,
  method = "rpart"
)

dt_importance <- varImp(
  dt_caret_model,
  scale = TRUE
)

dt_feature_df <- data.frame(
  Feature = rownames(dt_importance$importance),
  Importance = dt_importance$importance[, 1]
)

rownames(dt_feature_df) <- NULL

View(dt_feature_df)

ggplot(
  dt_feature_df,
  aes(
    x = reorder(Feature, Importance),
    y = Importance,
    fill = Importance
  )
) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "gold3",
    high = "maroon4"
  ) +
  coord_flip() +
  labs(
    title = "Feature Importance in Decision Tree",
    x = "Features",
    y = "Importance"
  )
