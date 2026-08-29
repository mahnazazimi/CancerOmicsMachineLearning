data <- read.table(
  "Raw Matrix.txt",
  sep = "\t"
)

# Remove genes with zero expression across all samples
data <- data[rowSums(data) > 0, ]

# CPM normalization
data <- apply(
  data,
  2,
  function(x) (x / sum(x)) * 1e6
)

# Differentially expressed genes
DEG <- read.table(
  "DEG.txt",
  sep = "\t"
)

DEG <- subset(
  DEG,
  padj < 0.05 & abs(log2FoldChange) > 5
)

# Keep only selected DEGs
data <- data[
  rownames(data) %in% rownames(DEG),
]

# Reshape to machine learning data structure:
# samples as rows and genes as columns
data <- t(data)
data <- as.data.frame(data)

# Metadata preparation
metadata <- read.table(
  "MetaData.txt",
  sep = "\t",
  header = TRUE
)

metadata$SampleName <- gsub(
  "-",
  ".",
  metadata$SampleName
)

# Confirm sample order between expression data and metadata
all(
  rownames(data) == metadata$SampleName
)


# Add classification labels
metadata$SampleType <- factor(
  metadata$SampleType,
  levels = c("Normal", "Tumor")
)

data$Label <- metadata$SampleType


# Stratified train/test split
library(caret)

train_index <- createDataPartition(
  data$Label,
  p = 0.8,
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

table(train_data$Label)
table(test_data$Label)


# Logistic Regression
lr_model <- glm(
  Label ~ .,
  data = train_data,
  family = binomial
)

lr_probability <- predict(
  lr_model,
  newdata = test_data,
  type = "response"
)

lr_pred <- ifelse(
  lr_probability >= 0.5,
  "Tumor",
  "Normal"
)

lr_pred <- factor(
  lr_pred,
  levels = levels(data$Label)
)

lr_result <- data.frame(
  Actual = test_data$Label,
  Predicted = lr_pred
)


# Support Vector Machine
library(e1071)

svm_model <- svm(
  Label ~ .,
  data = train_data
)

svm_pred <- predict(
  svm_model,
  newdata = test_data
)

svm_result <- data.frame(
  Actual = test_data$Label,
  Predicted = svm_pred
)


# Decision Tree
library(rpart)

dt_model <- rpart(
  Label ~ .,
  data = train_data,
  method = "class"
)

dt_pred <- predict(
  dt_model,
  newdata = test_data,
  type = "class"
)

dt_result <- data.frame(
  Actual = test_data$Label,
  Predicted = dt_pred
)


# Random Forest
library(randomForest)

rf_model <- randomForest(
  Label ~ .,
  data = train_data
)

rf_pred <- predict(
  rf_model,
  newdata = test_data
)

rf_result <- data.frame(
  Actual = test_data$Label,
  Predicted = rf_pred
)


# Naive Bayes
library(e1071)

nb_model <- naiveBayes(
  Label ~ .,
  data = train_data
)

nb_pred <- predict(
  nb_model,
  newdata = test_data
)

nb_result <- data.frame(
  Actual = test_data$Label,
  Predicted = nb_pred
)


# K-Nearest Neighbors
library(caret)

knn_model <- train(
  Label ~ .,
  data = train_data,
  method = "knn"
)

knn_pred <- predict(
  knn_model,
  newdata = test_data
)

knn_result <- data.frame(
  Actual = test_data$Label,
  Predicted = knn_pred
)


# Store trained models
all_models <- list(
  Logistic_Regression = lr_model,
  SVM = svm_model,
  Decision_Tree = dt_model,
  Random_Forest = rf_model,
  Naive_Bayes = nb_model,
  KNN = knn_model
)
