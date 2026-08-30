# STAD Candidate Biomarker Identification

library(caret)
library(randomForest)
library(ggplot2)

set.seed(123)


count_data <- read.table(
  "STAD/Raw Matrix.txt",
  sep = "\t"
)

deg_data <- read.table(
  "STAD/DEG.txt",
  sep = "\t"
)

metadata <- read.table(
  "STAD/MetaData.txt",
  sep = "\t",
  header = TRUE
)

# DEG filtering

deg_filtered <- subset(
  deg_data,
  padj < 0.05 & abs(log2FoldChange) > 2
)

# Keep expression values only for selected DEGs
count_data <- count_data[
  rownames(count_data) %in% rownames(deg_filtered),
]

# Reshape to machine learning data structure
# samples as rows, genes as columns

count_data <- t(count_data)
count_data <- as.data.frame(count_data)

metadata$SampleName <- gsub(
  "-",
  ".",
  metadata$SampleName
)

# Confirm that expression samples and metadata match
stopifnot(
  all(rownames(count_data) == metadata$SampleName)
)


# Convert expression features to numeric
count_data[] <- lapply(
  count_data,
  as.numeric
)

# Add tumor/normal labels
count_data$Labels <- factor(
  metadata$SampleType
)

data <- count_data

# Train/test split

train_index <- createDataPartition(
  data$Labels,
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
    Labels ~ .,
    data = train_data,
    method = "glm",
    trControl = train_control
  ),

  "SVM" = train(
    Labels ~ .,
    data = train_data,
    method = "svmRadial",
    trControl = train_control
  ),

  "Decision Tree" = train(
    Labels ~ .,
    data = train_data,
    method = "rpart",
    trControl = train_control
  ),

  "Random Forest" = train(
    Labels ~ .,
    data = train_data,
    method = "rf",
    trControl = train_control
  ),

  "Naive Bayes" = train(
    Labels ~ .,
    data = train_data,
    method = "naive_bayes",
    trControl = train_control
  ),

  "KNN" = train(
    Labels ~ .,
    data = train_data,
    method = "knn",
    trControl = train_control
  )
)

cv_results <- resamples(models)

summary(cv_results)

dotplot(
  cv_results,
  main = "STAD Model Comparison: 5-Fold Cross-Validation"
)

# Random Forest classification

rf_model <- randomForest(
  Labels ~ .,
  data = train_data,
  importance = TRUE,
  ntree = 100
)

rf_pred <- predict(
  rf_model,
  newdata = test_data
)

conf_matrix <- confusionMatrix(
  rf_pred,
  test_data$Labels
)

conf_matrix

# Confusion Matrix visualization

conf_matrix_data <- as.data.frame(
  as.table(conf_matrix$table)
)

ggplot(
  conf_matrix_data,
  aes(
    Prediction,
    Reference,
    fill = Freq
  )
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
    title = "Random Forest Confusion Matrix"
  ) +
  theme_minimal()


# Gene feature importance

importance_df <- as.data.frame(
  importance(rf_model)
)

importance_df$Feature <- rownames(
  importance_df
)

importance_df <- importance_df[
  order(
    importance_df$MeanDecreaseGini,
    decreasing = TRUE
  ),
]

rownames(importance_df) <- NULL


# Feature importance visualization

ggplot(
  importance_df,
  aes(
    x = reorder(Feature, MeanDecreaseGini),
    y = MeanDecreaseGini,
    fill = MeanDecreaseGini
  )
) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "gold3",
    high = "maroon4"
  ) +
  coord_flip() +
  labs(
    title = "Gene Importance in Random Forest",
    x = "Gene Features",
    y = "Mean Decrease Gini"
  )


# Map gene IDs to gene symbols

annotation <- read.table(
  "annot_df.txt",
  sep = "\t"
)

importance_df$SYMBOL <- annotation[
  match(
    importance_df$Feature,
    annotation$ID
  ),
  2
]


# Export ranked candidate biomarkers

write.table(
  importance_df,
  file = "stad_candidate_biomarkers.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
