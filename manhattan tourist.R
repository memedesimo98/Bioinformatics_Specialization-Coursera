# manhattam tourist

manhattan_tourist <- function(n, m, Down, Right) {
  s <- matrix(0, n+1, m+1)  # Initialize DP table
  
  for (i in 1:(n+1)) {
    for (j in 1:(m+1)) {
      if (i > 1) s[i, j] <- max(s[i, j], s[i-1, j] + Down[i-1, j])
      if (j > 1) s[i, j] <- max(s[i, j], s[i, j-1] + Right[i, j-1])
    }
  }
  
  return(s[n+1, m+1])  # Final result at (n,m)
}

# Example Input Data
n <- 11
m <- 16
Down <- matrix(c(
  4,0,3,1,0,2,2,3,1,4,1,2,1,2,0,2,2,
  3,3,1,2,0,2,3,4,0,1,1,1,4,2,0,0,4,
  4,2,3,0,1,3,0,4,0,2,0,4,1,4,1,0,1,
  4,0,1,0,4,2,0,4,2,0,2,4,4,2,2,1,4,
  0,0,2,0,1,3,0,4,0,4,2,4,1,0,4,0,2,
  3,1,3,0,1,0,2,3,4,0,3,3,3,0,4,2,1,
  2,0,4,1,0,0,1,1,1,1,4,1,4,1,0,0,0,
  2,4,2,1,1,3,0,4,2,2,1,3,1,0,4,0,3,
  4,4,0,0,4,0,0,4,0,0,1,2,2,3,4,2,4,
  1,0,3,4,4,2,3,1,0,3,4,0,1,4,3,2,1,
  4,4,4,0,2,4,4,0,4,3,0,3,2,4,4,4,2
), n, m+1, byrow=TRUE)
Right <- matrix(c(
  4,2,4,4,4,3,1,4,2,4,2,0,1,0,3,4,
  4,1,0,0,4,4,3,2,3,3,2,2,3,0,4,3,
  0,4,1,1,2,4,0,4,3,4,4,3,0,0,1,2,
  0,0,0,0,2,1,1,0,3,1,0,4,3,2,2,4,
  0,3,0,1,0,3,4,0,1,1,0,3,4,2,4,4,
  0,0,4,1,2,3,1,1,0,3,0,4,1,1,1,2,
  2,2,0,4,0,1,4,3,1,4,2,1,0,0,3,4,
  4,4,0,4,2,4,3,3,2,4,4,0,4,3,0,4,
  3,4,0,0,0,4,3,4,1,0,0,1,3,0,4,2,
  1,3,2,2,0,2,2,3,4,3,3,4,3,2,3,0,
  3,3,2,3,2,2,0,0,4,2,3,0,1,3,4,1,
  0,4,2,0,0,4,2,3,3,3,1,3,1,1,0,1
), n+1, m, byrow=TRUE)

manhattan_tourist(n, m, Down, Right)  

v <- "ABBA"


LCS_backtracer <- function(v, w) {
  v <- unlist(strsplit(v, ""))
  w <- unlist(strsplit(w, ""))
  V <- length(v)
  W <- length(w)
  
  s <- matrix(0, nrow = V+1, ncol = W+1)
  Backtrack <- matrix("", nrow = V+1, ncol = W+1)
  
  for (i in 2:(V+1)) {
    for (j in 2:(W+1)) {
      match <- ifelse(v[i-1] == w[j-1], 1, 0)
      s[i, j] <- max(s[i-1, j], s[i, j-1], s[i-1, j-1] + match)
      
      if (s[i, j] == s[i-1, j-1] + match) {
        Backtrack[i, j] <- "dr"
      } else if (s[i, j] == s[i-1, j]) {
        Backtrack[i, j] <- "d"
      } else {
        Backtrack[i, j] <- "r"
      }
    }
  }
  
  alignment <- Output_LCS(Backtrack,v,V+1,W+1)
  return(list(Backtrack = Backtrack, matrix = s,End_score = s[V+1,W+1], Alignment = alignment, 
              Edit_distance = (V + W) - 2 * s[V+1, W+1]) )
}

Output_LCS <- function(Backtrack, v, i, j) {
  if (i == 0 || j == 0) {
    return("")
  }
  if (Backtrack[i, j] == "d") {
    return(Output_LCS(Backtrack, v, i-1, j))
  }
  if (Backtrack[i, j] == "r") {
    return(Output_LCS(Backtrack, v, i, j-1))
  }
  if (Backtrack[i, j] == "dr") {
    return(paste0(Output_LCS(Backtrack, v, i-1, j-1), v[i-1]))
  }
}

v <- "TGTACG"
w <- "GCTAGT"
S <- LCS_backtracer(v,w)

vv <- unlist(strsplit(v, ""))
ww <- unlist(strsplit(w, ""))
V <- length(vv)
W <- length(ww)

SS <- Output_LCS(S,vv,V,W)
