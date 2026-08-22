# Regularization on TCGA-ESCA gene-expression data

library(caret)
library(glmnet)

data <- read.table("data.txt", sep = "\t", header = TRUE)
data <- t(data)  # Reshape to machine learning data structure: samples as rows, genes as columns
data <- as.data.frame(data)

train_index <- createDataPartition(
  data$G6PD,
  p = 0.8,
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

response_gene <- "G6PD"
predictor_genes <- colnames(train_data)[2:10]

x <- as.matrix(train_data[, predictor_genes])
y <- train_data[[response_gene]]

# Standard multiple linear regression
model_formula <- as.formula(
  paste(response_gene, "~", paste(predictor_genes, collapse = " + "))
)

linear_model <- lm(
  model_formula,
  data = train_data
)

coef(linear_model)


# Ridge regression
# alpha = 0 -> Ridge
# alpha = 1 -> Lasso
# 0 < alpha < 1 -> Elastic Net

ridge_model <- glmnet(
  x,
  y,
  alpha = 0
)

# Inspect coefficients at a selected lambda
coef(
  ridge_model,
  s = 0.1
)

# Cross-validation to select the lambda with minimum CV error
cv_ridge <- cv.glmnet(
  x,
  y,
  alpha = 0
)

plot(cv_ridge)

best_lambda_ridge <- cv_ridge$lambda.min

coef(
  cv_ridge,
  s = best_lambda_ridge
)


# Lasso regression for feature selection
lasso_model <- glmnet(
  x,
  y,
  alpha = 1
)

# Some coefficients may shrink to zero
coef(
  lasso_model,
  s = 0.1
)

# Cross-validation to select the optimal lambda
cv_lasso <- cv.glmnet(
  x,
  y,
  alpha = 1
)

plot(cv_lasso)

best_lambda_lasso <- cv_lasso$lambda.min

lasso_coefficients <- coef(
  cv_lasso,
  s = best_lambda_lasso
)

lasso_coefficients
