# Breast Cancer Clustering with K-means

library(Rtsne)
library(ggplot2)

set.seed(123)

data <- read.csv("breast-cancer.csv")

rownames(data) <- data$id

diagnosis <- data$diagnosis

# Keep numerical diagnostic features only
features <- data[, !colnames(data) %in% c("id", "diagnosis")]

# Standardize features before distance-based clustering
features_scaled <- scale(features)

# Find an appropriate number of clusters using Elbow Method

max_k <- 10
wss <- numeric(max_k)

for (k in 1:max_k) {

  km_result <- kmeans(
    features_scaled,
    centers = k,
    nstart = 25
  )

  wss[k] <- km_result$tot.withinss
}

elbow_data <- data.frame(
  k = 1:max_k,
  WSS = wss
)

ggplot(
  elbow_data,
  aes(x = k, y = WSS)
) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 1:max_k) +
  labs(
    title = "Elbow Plot for K-means Clustering",
    x = "Number of Clusters (k)",
    y = "Within-Cluster Sum of Squares"
  ) +
  theme_minimal()

# K-means clustering

kmeans_result <- kmeans(
  features_scaled,
  centers = 2,
  nstart = 25
)

cluster_labels <- kmeans_result$cluster

# t-SNE visualization

tsne_result <- Rtsne(
  features_scaled,
  dims = 2
)

tsne_data <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  Cluster = as.factor(cluster_labels),
  Diagnosis = diagnosis
)


# Visualize K-means clusters

ggplot(
  tsne_data,
  aes(
    x = tSNE1,
    y = tSNE2,
    color = Cluster
  )
) +
  geom_point() +
  labs(
    title = "K-means Clustering of Breast Cancer Data",
    x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_minimal()


# Visualize known diagnostic labels for comparison

ggplot(
  tsne_data,
  aes(
    x = tSNE1,
    y = tSNE2,
    color = Diagnosis
  )
) +
  geom_point() +
  labs(
    title = "Breast Cancer Diagnostic Labels",
    x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  theme_minimal()
