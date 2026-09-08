# Advanced Clustering Methods
# on Breast Cancer Diagnostic Data

library(kernlab)
library(mclust)
library(Rtsne)

set.seed(123)

data <- read.csv("breast-cancer.csv")

rownames(data) <- data$id
diagnosis <- data$diagnosis

features <- data[, !colnames(data) %in% c("id", "diagnosis")]

features_scaled <- scale(features)

# t-SNE representation for visualization

tsne_result <- Rtsne(
  features_scaled,
  dims = 2
)

tsne_data <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2]
)

# Spectral Clustering

spectral_model <- specc(
  as.matrix(features_scaled),
  centers = 2
)

spectral_clusters <- as.integer(
  spectral_model
)

plot(
  tsne_data$tSNE1,
  tsne_data$tSNE2,
  col = spectral_clusters,
  pch = 19,
  xlab = "t-SNE 1",
  ylab = "t-SNE 2",
  main = "Spectral Clustering"
)

# Gaussian Mixture Model

gmm_model <- Mclust(
  features_scaled,
  G = 2
)

gmm_clusters <- gmm_model$classification

summary(gmm_model)

plot(
  tsne_data$tSNE1,
  tsne_data$tSNE2,
  col = gmm_clusters,
  pch = 19,
  xlab = "t-SNE 1",
  ylab = "t-SNE 2",
  main = "Gaussian Mixture Model Clustering"
)
