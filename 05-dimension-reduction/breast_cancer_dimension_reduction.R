# Dimension Reduction on Breast Cancer Diagnostic Data

library(uwot)
library(Rtsne)

set.seed(123)

data <- read.csv("breast-cancer.csv")

rownames(data) <- data$id

labels <- data$diagnosis

# Remove ID and diagnosis columns
features <- data[, !colnames(data) %in% c("id", "diagnosis")]

# Standardize numerical features
features_scaled <- scale(features)

# Colors for visualization
sample_colors <- ifelse(
  labels == "M",
  "tomato",
  "seagreen"
)

# Principal Component Analysis (PCA)

pca_model <- prcomp(
  features_scaled,
  center = TRUE
)

pca_data <- data.frame(
  PC1 = pca_model$x[, 1],
  PC2 = pca_model$x[, 2],
  Diagnosis = labels
)

plot(
  pca_data$PC1,
  pca_data$PC2,
  col = sample_colors,
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA"
)

# UMAP

umap_result <- umap(
  features_scaled
)

umap_data <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Diagnosis = labels
)

plot(
  umap_data$UMAP1,
  umap_data$UMAP2,
  col = sample_colors,
  pch = 19,
  cex = 0.8,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = "UMAP"
)

# t-SNE

tsne_result <- Rtsne(
  features_scaled
)

tsne_data <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  Diagnosis = labels
)

plot(
  tsne_data$tSNE1,
  tsne_data$tSNE2,
  col = sample_colors,
  pch = 19,
  cex = 0.8,
  xlab = "t-SNE1",
  ylab = "t-SNE2",
  main = "t-SNE"
)
