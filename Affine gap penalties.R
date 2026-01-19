
### Affine gap penalties

Affine_Global_Alignment <- function(v, w, match = 1, mismatch = 1, gap_open = 2, gap_extend = 1) {
  v <- unlist(strsplit(v, "")) 
  w <- unlist(strsplit(w, "")) 
  V <- length(v)
  W <- length(w)
  
  # Initialize score matrices
  M <- matrix(-Inf, nrow = V+1, ncol = W+1)
  I <- matrix(-Inf, nrow = V+1, ncol = W+1)
  D <- matrix(-Inf, nrow = V+1, ncol = W+1)
  Backtrack <- matrix("", nrow = V+1, ncol = W+1)
  
  # Initialization
  M[1,1] <- 0
  for (i in 2:(V+1)) {
    D[i,1] <- -gap_open - (i - 2) * gap_extend
    M[i,1] <- D[i,1]
    Backtrack[i,1] <- "d"
  }
  for (j in 2:(W+1)) {
    I[1,j] <- -gap_open - (j - 2) * gap_extend
    M[1,j] <- I[1,j]
    Backtrack[1,j] <- "r"
  }
  
  # Fill matrices
  for (i in 2:(V+1)) {
    for (j in 2:(W+1)) {
      # Match/mismatch
      score_diag <- ifelse(v[i-1] == w[j-1], match, -mismatch)
      from_M <- M[i-1, j-1] + score_diag
      from_I <- I[i-1, j-1] + score_diag
      from_D <- D[i-1, j-1] + score_diag
      M[i,j] <- max(from_M, from_I, from_D)
      
      # Insertion: gap in v (move right)
      I[i,j] <- max(
        M[i, j-1] - gap_open,
        I[i, j-1] - gap_extend
      )
      
      # Deletion: gap in w (move down)
      D[i,j] <- max(
        M[i-1, j] - gap_open,
        D[i-1, j] - gap_extend
      )
      
      # Choose best overall path
      scores <- c(M[i,j], I[i,j], D[i,j])
      direction <- c("dr", "r", "d")
      best <- which.max(scores)
      Backtrack[i, j] <- direction[best]
    }
  }
  
  # Traceback
  Output_Alignment <- function(Backtrack, v, w, i, j) {
    if (i == 1 && j == 1) return(list(v_align = "", w_align = ""))
    
    if (i > 1 && Backtrack[i, j] == "d") {
      res <- Output_Alignment(Backtrack, v, w, i - 1, j)
      return(list(v_align = paste0(res$v_align, v[i - 1]),
                  w_align = paste0(res$w_align, "-")))
    }
    
    if (j > 1 && Backtrack[i, j] == "r") {
      res <- Output_Alignment(Backtrack, v, w, i, j - 1)
      return(list(v_align = paste0(res$v_align, "-"),
                  w_align = paste0(res$w_align, w[j - 1])))
    }
    
    if (i > 1 && j > 1 && Backtrack[i, j] == "dr") {
      res <- Output_Alignment(Backtrack, v, w, i - 1, j - 1)
      return(list(v_align = paste0(res$v_align, v[i - 1]),
                  w_align = paste0(res$w_align, w[j - 1])))
    }
    
    return(list(v_align = "", w_align = ""))
  }
  
  alignment <- Output_Alignment(Backtrack, v, w, V+1, W+1)
  final_score <- max(M[V+1, W+1], I[V+1, W+1], D[V+1, W+1])
  
  return(list(
    Score = final_score,
    Alignment = alignment,
    Backtrack = Backtrack
  ))
}

v <- "GA"
w <-"GTTA"
result <- Affine_Global_Alignment(v,w,1,3,2,1)

Affine_Local_Alignment <- function(v, w, match = 1, mismatch = 1, gap_open = 2, gap_extend = 1) {
  v <- unlist(strsplit(v, ""))
  w <- unlist(strsplit(w, ""))
  V <- length(v)
  W <- length(w)
  
  # Initialize matrices
  M <- matrix(0, nrow = V + 1, ncol = W + 1)
  I <- matrix(0, nrow = V + 1, ncol = W + 1)
  D <- matrix(0, nrow = V + 1, ncol = W + 1)
  Backtrack <- matrix("", nrow = V + 1, ncol = W + 1)
  
  # Track max score and its position
  max_score <- 0
  max_i <- 1
  max_j <- 1
  
  for (i in 2:(V + 1)) {
    for (j in 2:(W + 1)) {
      # Diagonal: match or mismatch
      score <- ifelse(v[i - 1] == w[j - 1], match, -mismatch)
      from_M <- M[i - 1, j - 1] + score
      from_I <- I[i - 1, j - 1] + score
      from_D <- D[i - 1, j - 1] + score
      M[i, j] <- max(0, from_M, from_I, from_D)
      
      # Right: gap in v
      I[i, j] <- max(
        M[i, j - 1] - gap_open,
        I[i, j - 1] - gap_extend,
        0
      )
      
      # Down: gap in w
      D[i, j] <- max(
        M[i - 1, j] - gap_open,
        D[i - 1, j] - gap_extend,
        0
      )
      
      # Backtrack with priority: dr > r > d
      scores <- c(M[i, j], I[i, j], D[i, j])
      directions <- c("dr", "r", "d")
      best <- which.max(scores)
      Backtrack[i, j] <- directions[best]
      
      # Update max score
      if (scores[best] > max_score) {
        max_score <- scores[best]
        max_i <- i
        max_j <- j
      }
    }
  }
  
  # Traceback until score hits 0
  Output_Alignment <- function(Backtrack, M, v, w, i, j) {
    v_align <- ""
    w_align <- ""
    while (i > 1 || j > 1) {
      if (M[i, j] == 0) break
      if (Backtrack[i, j] == "dr") {
        v_align <- paste0(v[i - 1], v_align)
        w_align <- paste0(w[j - 1], w_align)
        i <- i - 1
        j <- j - 1
      } else if (Backtrack[i, j] == "r") {
        v_align <- paste0("-", v_align)
        w_align <- paste0(w[j - 1], w_align)
        j <- j - 1
      } else if (Backtrack[i, j] == "d") {
        v_align <- paste0(v[i - 1], v_align)
        w_align <- paste0("-", w_align)
        i <- i - 1
      } else {
        break
      }
    }
    return(list(v_align = v_align, w_align = w_align))
  }
  
  alignment <- Output_Alignment(Backtrack, M, v, w, max_i, max_j)
  return(list(
    Score = max_score,
    Alignment = alignment,
    AlignmentEndsAt = list(i = max_i, j = max_j)
  ))
}
