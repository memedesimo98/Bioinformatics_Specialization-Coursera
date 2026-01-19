LinearSpaceAlignment <- function(v, w, top, bottom, left, right, match = 1, mismatch = 1, indel = 1) {
  # Base cases
  if (left == right) {
    # Vertical edges: align v[top+1 ... bottom] to empty string w
    return(rep("V", bottom - top))
  }
  if (top == bottom) {
    # Horizontal edges: align w[left+1 ... right] to empty string v
    return(rep("H", right - left))
  }
  
  # Compute middle edge
  middle <- floor((left + right) / 2)
  mid_info <- Find_Middle_Edge(
    substr(v, top + 1, bottom),
    substr(w, left + 1, right),
    match, mismatch, indel
  )
  
  # Extract middle edge positions
  i_rel <- mid_info$from[1]
  j_rel <- mid_info$from[2]
  i_abs <- i_rel + top
  j_abs <- j_rel + left
  
  move <- mid_info$move
  mid_move <- switch(move, "↓" = "V", "→" = "H", "↘" = "D")
  
  # Adjust positions based on edge direction
  next_i <- if (move %in% c("↓", "↘")) i_abs + 1 else i_abs
  next_j <- if (move %in% c("→", "↘")) j_abs + 1 else j_abs
  
  # Recursively compute path
  left_path <- LinearSpaceAlignment(v, w, top, i_abs, left, j_abs, match, mismatch, indel)
  right_path <- LinearSpaceAlignment(v, w, next_i, bottom, next_j, right, match, mismatch, indel)
  
  return(c(left_path, mid_move, right_path))
}

BuildAlignment <- function(v, w, path, match = 1, mismatch = 1, indel = 1) {
  v <- unlist(strsplit(v, ""))
  w <- unlist(strsplit(w, ""))
  i <- j <- 1
  v_aligned <- c()
  w_aligned <- c()
  score <- 0
  
  for (move in path) {
    if (move == "D") {
      v_aligned <- c(v_aligned, v[i])
      w_aligned <- c(w_aligned, w[j])
      if (v[i] == w[j]) {
        score <- score + match
      } else {
        score <- score - mismatch
      }
      i <- i + 1
      j <- j + 1
    } else if (move == "H") {
      v_aligned <- c(v_aligned, "-")
      w_aligned <- c(w_aligned, w[j])
      score <- score - indel
      j <- j + 1
    } else if (move == "V") {
      v_aligned <- c(v_aligned, v[i])
      w_aligned <- c(w_aligned, "-")
      score <- score - indel
      i <- i + 1
    }
  }
  
  return(list(
    End_score = score,
    Alignment = list(
      v_align = paste(v_aligned, collapse = ""),
      w_align = paste(w_aligned, collapse = "")
    )
  ))
}


# Step 1: generate alignment path
path <- LinearSpaceAlignment(v, w, top = 0, bottom = nchar(v), left = 0, right = nchar(w), match = 1, mismatch = 1, indel = 5)

# Step 2: build full alignment and score
result_2 <- BuildAlignment(v, w, path, match = 1, mismatch = 1, indel = 5)

# Output
print(result_2$End_score)
cat(result_2$Alignment$v_align, "\n")
cat(result_2$Alignment$w_align, "\n")
