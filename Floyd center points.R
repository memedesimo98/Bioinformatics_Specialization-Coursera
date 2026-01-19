### Floyd center points 

# Euclidean distance function
euclidean <- function(a, b) sqrt(sum((a - b)^2))

Floyd_centers <- function(data,k,old_centers_ids = NULL) {
  n <- nrow(data)
  if (is.null(old_centers_ids)) {
  old_centers_ids <- sample(1:n,k)
  old_centers <- data[old_centers_ids,]
  } else {
    old_centers <- old_centers_ids
  } 
  new_center <- c()
  while (TRUE) {
    # save previous centers
    prev_centers <- old_centers
    
    # assignment step
    assignments <- vector("list", nrow(prev_centers))
    for (i in 1:nrow(data)) {
      dists <- sapply(1:nrow(prev_centers), function(j) euclidean(data[i, ], prev_centers[j, ]))
      nearest_center <- which.min(dists)
      assignments[[nearest_center]] <- rbind(assignments[[nearest_center]], data[i, ])
    }
    
    # update step
    gravity_points <- lapply(assignments, function(cluster) {
      if (is.null(cluster)) return(NULL) else colMeans(cluster)
    })
    new_center <- do.call(rbind, gravity_points)
    
    # Round to 3 decimals for both reporting and convergence check
    new_center <- round(new_center, 3)
    
    if (identical(new_center, old_centers)) {
      break
    }
    
    old_centers <- new_center
  }
  
  return(new_center)
}

old_centers_ids <- c(1:k)
old_centers <- data[old_centers_ids,]
centers<-Floyd_centers(data_matrix,k,old_centers_ids)

apply(centers, 1, function(row) {
  cat(paste(row, collapse = " "), "\n")
})