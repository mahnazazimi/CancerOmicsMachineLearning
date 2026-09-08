# Density-Based and Hierarchical Clustering
# on Breast Cancer Diagnostic Data

library(dbscan)
library(Rtsne)

set.seed(123)

# Data preparation

data <- read.csv("breast-cancer.csv")

rownames(data) <- data$id
diagnosis <- data$diagnosis

features <- data[, !colnames(data) %in% c("id", "diagnosis")]

features_scaled <- scale(features)

# DBSCAN

# Use t-SNE representation for exploratory DBSCAN visualization

tsne_result <- Rtsne(
  features_scaled,
  dims = 2
)

tsne_data <- data.frame(
  Dim1 = tsne_result$Y[, 1],
  Dim2 = tsne_result$Y[, 2]
)

# Inspect nearest-neighbor distances to help choose epsilon

k <- 4

kNNdistplot(
  tsne_data,
  k = k
)

# Set epsilon after inspecting the kNN distance plot
eps_value <- 1.6

dbscan_result <- dbscan(
  tsne_data,
  eps = eps_value,
  minPts = 4
)

dbscan_labels <- dbscan_result$cluster

table(dbscan_labels)

plot(
  tsne_data$Dim1,
  tsne_data$Dim2,
  col = dbscan_labels + 1,
  pch = 19,
  xlab = "t-SNE 1",
  ylab = "t-SNE 2",
  main = "DBSCAN Clustering"
)

# Hierarchical Clustering

distance_matrix <- dist(
  features_scaled
)

hc_model <- hclust(
  distance_matrix,
  method = "ward.D2"
)

plot(
  hc_model,
  labels = FALSE,
  main = "Hierarchical Clustering Dendrogram",
  xlab = "Samples",
  ylab = "Distance"
)

hierarchical_clusters <- cutree(
  hc_model,
  k = 2
)

table(hierarchical_clusters)
