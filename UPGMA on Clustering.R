### UPGMA on Clustering

matrix_to_clusters <- function(D, labels = NULL, method = "average") {
  # Step 1: Perform hierarchical clustering
  hc <- hclust(as.dist(D), method = method)
  
  n <- length(hc$order)
  clusters <- vector("list", n - 1)
  
  # Step 2: Reconstruct clusters from hc$merge
  for (i in seq_len(n - 1)) {
    members <- c()
    for (j in 1:2) {
      idx <- hc$merge[i, j]
      if (idx < 0) {
        # leaf
        members <- c(members, -idx)
      } else {
        # previously formed cluster
        members <- c(members, clusters[[idx]])
      }
    }
    clusters[[i]] <- sort(members)
  }
  
  # Step 3: Apply labels if provided
  if (!is.null(labels)) {
    clusters <- lapply(clusters, function(cl) labels[cl])
  }
  
  return(clusters)
}

D <- matrix(c(
  0.00, 0.74, 0.85, 0.54, 0.83, 0.92, 0.89,
  0.74, 0.00, 1.59, 1.35, 1.20, 1.48, 1.55,
  0.85, 1.59, 0.00, 0.63, 1.13, 0.69, 0.73,
  0.54, 1.35, 0.63, 0.00, 0.66, 0.43, 0.88,
  0.83, 1.20, 1.13, 0.66, 0.00, 0.72, 0.55,
  0.92, 1.48, 0.69, 0.43, 0.72, 0.00, 0.80,
  0.89, 1.55, 0.73, 0.88, 0.55, 0.80, 0.00
), nrow = 7, byrow = TRUE)

clusters <- matrix_to_clusters(D)
# Suppose clusters is your list of integer vectors
cat(sapply(clusters, function(x) paste(x, collapse = " ")), sep = "\n")

