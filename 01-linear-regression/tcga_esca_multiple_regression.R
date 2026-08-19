# Data preparation
data <- read.table("data.txt", sep = "\t", header = TRUE)
data <- t(data)  
data <- as.data.frame(data)

library(caret)
library(Metrics)

# Train/test split
train_index <- createDataPartition(
  data$G6PD,
  p = 0.8,
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]


# Multiple linear regression using the first 10 predictor genes
response_gene <- colnames(data)[1]
predictor_genes <- colnames(data)[2:11]

model_formula <- as.formula(
  paste(response_gene, "~", paste(predictor_genes, collapse = " + "))
)

model_extended <- lm(
  model_formula,
  data = train_data
)

predicted_extended <- predict(
  model_extended,
  newdata = test_data
)

actual <- test_data[[response_gene]]

rmse(actual, predicted_extended)
