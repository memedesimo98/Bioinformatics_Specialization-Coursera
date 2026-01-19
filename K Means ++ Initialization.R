### K Means ++ Initialization

# k-Means++ initializer
kmeanspp_init <- function(data, k) {
  n <- nrow(data)
  centers <- list()
  
  # 1. pick first center uniformly at random
  first_idx <- sample(1:n, 1)
  centers[[1]] <- data[first_idx, ]
  
  # 2. pick remaining centers
  while (length(centers) < k) {
    # compute squared distance of each point to nearest center
    dists <- sapply(1:n, function(i) {
      min(sapply(centers, function(c) euclidean(data[i, ], c)))
    })
    dists2 <- dists^2
    
    # normalize to probabilities
    probs <- dists2 / sum(dists2)
    
    # sample next center with probability proportional to squared distance
    next_idx <- sample(1:n, 1, prob = probs)
    centers[[length(centers) + 1]] <- data[next_idx, ]
  }
  
  do.call(rbind, centers)
}

# Euclidean distance function
euclidean <- function(a, b) sqrt(sum((a - b)^2))

init <- kmeanspp_init(data, k)
final <- Floyd_centers(data, k, init)
apply(final, 1, function(row) {
  cat(paste(row, collapse = " "), "\n")
})
