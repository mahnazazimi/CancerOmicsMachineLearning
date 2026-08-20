data <- read.table("data.txt", sep = "\t", header = TRUE)
data <- t(data)  # Reshape to machine learning data structure: samples as rows, genes as columns
data <- as.data.frame(data)

# Detailed matrix plot for the first 10 genes
library(GGally)

png(
  "matrix_plot.png",
  width = 4000,
  height = 4000,
  res = 300
)

ggpairs(data[, 1:10])
dev.off()


# Multiple regression with two predictors
data <- data[, 1:3]

colnames(data) <- c(
  "response",
  "predictor_1",
  "predictor_2"
)

# Train/test split
library(caret)

train_index <- createDataPartition(
  data$predictor_1,
  p = 0.8,
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

# Fit multiple linear regression using training data
model <- lm(
  response ~ predictor_1 + predictor_2,
  data = train_data
)


# Create a grid for the regression surface
grid_points <- expand.grid(
  predictor_1 = seq(
    min(c(train_data$predictor_1, test_data$predictor_1)),
    max(c(train_data$predictor_1, test_data$predictor_1)),
    length.out = 30
  ),
  
  predictor_2 = seq(
    min(c(train_data$predictor_2, test_data$predictor_2)),
    max(c(train_data$predictor_2, test_data$predictor_2)),
    length.out = 30
  )
)

# Predict response values across the grid
grid_points$response_pred <- predict(
  model,
  newdata = grid_points
)

# Predicted values for training and test samples
train_data$response_pred <- predict(
  model,
  newdata = train_data
)

test_data$response_pred <- predict(
  model,
  newdata = test_data
)


# Interactive 3D visualization
library(plotly)

p <- plot_ly() %>%
  
  add_markers(
    data = train_data,
    x = ~predictor_1,
    y = ~predictor_2,
    z = ~response,
    marker = list(color = "blue", size = 5),
    name = "Train"
  ) %>%
  
  add_markers(
    data = test_data,
    x = ~predictor_1,
    y = ~predictor_2,
    z = ~response,
    marker = list(color = "red", size = 5),
    name = "Test"
  ) %>%
  
  # Multiple regression is represented by a surface
  add_surface(
    x = matrix(grid_points$predictor_1, nrow = 30),
    y = matrix(grid_points$predictor_2, nrow = 30),
    z = matrix(grid_points$response_pred, nrow = 30),
    opacity = 0.5,
    showscale = FALSE,
    name = "Model"
  ) %>%
  
  layout(
    scene = list(
      title = "Multiple Linear Regression",
      xaxis = list(title = "Predictor 1"),
      yaxis = list(title = "Predictor 2"),
      zaxis = list(title = "Response")
    )
  )

p

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

# Define response and predictor features
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
