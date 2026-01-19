# middle edge
Find_Middle_Edge <- function(v, w, match = 1, mismatch = 1, indel = 1) {
  # Convert strings to character vectors
  v <- unlist(strsplit(v, ""))
  w <- unlist(strsplit(w, ""))
  
  
  if (length(v) > length(w)) {
    temp <- v
    v <- w
    w <- temp
  }
  
  n <- length(v)
  m <- length(w)
  mid <- floor(m / 2)
  
  score_fn <- function(a, b) ifelse(a == b, match, -mismatch)
  
  # Column-based forward DP: compute FromSource(i)
  compute_column_dp <- function(v, w_sub) {
    prev <- 0:n * -indel
    for (char in w_sub) {
      curr <- numeric(n + 1)
      curr[1] <- prev[1] - indel
      for (i in 2:(n + 1)) {
        diag <- prev[i - 1] + score_fn(v[i - 1], char)
        up   <- prev[i] - indel
        left <- curr[i - 1] - indel
        curr[i] <- max(diag, up, left)
      }
      prev <- curr
    }
    return(prev)
  }
  
  # Reverse DP: compute ToSink(i)
  compute_reverse_column_dp <- function(v, w_sub) {
    v_rev <- rev(v)
    w_rev <- rev(w_sub)
    prev <- 0:n * -indel
    for (char in w_rev) {
      curr <- numeric(n + 1)
      curr[1] <- prev[1] - indel
      for (i in 2:(n + 1)) {
        diag <- prev[i - 1] + score_fn(v_rev[i - 1], char)
        up   <- prev[i] - indel
        left <- curr[i - 1] - indel
        curr[i] <- max(diag, up, left)
      }
      prev <- curr
    }
    return(rev(prev))  # Align direction with forward DP
  }
  
  # Score vectors
  scoreL <- compute_column_dp(v, w[1:mid])
  scoreR <- compute_reverse_column_dp(v, w[(mid + 1):m])
  total <- scoreL + scoreR
  
  i <- which.max(total) - 1  # best row in middle column (0-based)
  
  # Determine the best move from (i, mid)
  diag <- if (i < n && mid < m)
    score_fn(v[i + 1], w[mid + 1]) else -Inf
  down <- -indel
  right <- -indel
  
  diag_score <- scoreL[i + 1] + diag + scoreR[i + 2]
  down_score <- scoreL[i + 2] + down + scoreR[i + 2]
  right_score <- scoreL[i + 1] + right + scoreR[i + 1]
  
  move <- which.max(c(diag_score, down_score, right_score))
  
  edge_moves <- list(
    c(i, mid, i + 1, mid + 1),  # Diagonal
    c(i, mid, i + 1, mid),      # Down
    c(i, mid, i, mid + 1)       # Right
  )
  
  return(list(
    from = edge_moves[[move]][1:2],
    to   = edge_moves[[move]][3:4],
    move = c("↘", "↓", "→")[move]
  ))
}





Find_Middle_Edge(v, w,1,1,2)
v <- "GAGA"
w <- "GAT"
