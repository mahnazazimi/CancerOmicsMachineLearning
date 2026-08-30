# Breast Cancer Feature Reduction

library(caret)
library(ggplot2)

set.seed(123)

data <- read.csv("breast-cancer.csv")

# Use sample IDs as row names and remove the ID column
rownames(data) <- data[, 1]
data <- data[, -1]

# Convert diagnosis to a classification label
data$diagnosis <- factor(
  data$diagnosis,
  levels = c("B", "M"),
  labels = c("Benign", "Malignant")
)


# Train/test split

train_index <- createDataPartition(
  data$diagnosis,
  p = 0.8,
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

# Compare classification algorithms using 5-fold CV

train_control <- trainControl(
  method = "cv",
  number = 5
)

models <- list(
  "Logistic Regression" = train(
    diagnosis ~ .,
    data = train_data,
    method = "glm",
    trControl = train_control
  ),

  "SVM" = train(
    diagnosis ~ .,
    data = train_data,
    method = "svmRadial",
    trControl = train_control
  ),

  "Decision Tree" = train(
    diagnosis ~ .,
    data = train_data,
    method = "rpart",
    trControl = train_control
  ),

  "Random Forest" = train(
    diagnosis ~ .,
    data = train_data,
    method = "rf",
    trControl = train_control
  ),

  "Naive Bayes" = train(
    diagnosis ~ .,
    data = train_data,
    method = "naive_bayes",
    trControl = train_control
  ),

  "KNN" = train(
    diagnosis ~ .,
    data = train_data,
    method = "knn",
    trControl = train_control
  )
)

cv_results <- resamples(models)

summary(cv_results)

dotplot(
  cv_results,
  main = "5-Fold Cross-Validation"
)

# SVM model using all features

svm_model_all <- train(
  diagnosis ~ .,
  data = train_data,
  method = "svmRadial"
)

pred_all <- predict(
  svm_model_all,
  newdata = test_data
)

conf_matrix_all <- confusionMatrix(
  pred_all,
  test_data$diagnosis,
  positive = "Malignant"
)

conf_matrix_all


# Feature importance

feature_importance <- varImp(
  svm_model_all,
  scale = TRUE
)

feature_df <- data.frame(
  Feature = rownames(feature_importance$importance),
  Importance = feature_importance$importance[, 1]
)

feature_df <- feature_df[
  order(feature_df$Importance, decreasing = TRUE),
]

rownames(feature_df) <- NULL


# Plot feature importance
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
    title = "Feature Importance in SVM",
    x = "Features",
    y = "Importance"
  )

# Feature reduction: select top 10 features

top_features <- feature_df$Feature[1:10]

reduced_train <- train_data[
  , c(top_features, "diagnosis")
]

reduced_test <- test_data[
  , c(top_features, "diagnosis")
]

# Train SVM using reduced feature set

svm_model_reduced <- train(
  diagnosis ~ .,
  data = reduced_train,
  method = "svmRadial"
)

pred_reduced <- predict(
  svm_model_reduced,
  newdata = reduced_test
)

conf_matrix_reduced <- confusionMatrix(
  pred_reduced,
  reduced_test$diagnosis,
  positive = "Malignant"
)

conf_matrix_reduced

# Compare full and reduced models

conf_matrix_all
conf_matrix_reduced
