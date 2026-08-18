data <- read.table("data.txt", sep = "\t")
data <- t(data)                    # Transpose to machine learning format: samples as rows, genes as columns
data <- as.data.frame(data)

predictor_gene <- "GCLC"
response_gene  <- "GCLM"

# Train/test split
library(caret)

train_index <- createDataPartition(
  data[[response_gene]],
  p = 0.8,                         # Use approximately 80% of samples for training
  list = FALSE
)

train_data <- data[train_index, ]
test_data  <- data[-train_index, ] # Remaining samples are kept for testing

# Visualize training and test data
plot(
  train_data[[predictor_gene]],
  train_data[[response_gene]],
  col = "seagreen",
  pch = 16,
  xlab = predictor_gene,
  ylab = response_gene
)

points(
  test_data[[predictor_gene]],
  test_data[[response_gene]],
  col = "tomato",
  pch = 16
)

# Train linear regression model
model_formula <- as.formula(
  paste(response_gene, "~", predictor_gene)
)                                 # Build formula dynamically from selected genes

model <- lm(
  y ~ x,
  data = train_data                # Fit the model using training data 
)

abline(model, col = "purple4")     # Add the fitted regression line to the plot

# Model coefficients
intercept <- coef(model)[1]        # Intercept of the fitted regression line
slope <- coef(model)[2]            # Effect of predictor on response

# Prediction on held-out test data
pred <- predict(
  model,
  newdata = test_data
)

# Compare predicted and observed expression values
result <- data.frame(
  predicted = pred,
  actual = test_data[[response_gene]]
)

View(result)
