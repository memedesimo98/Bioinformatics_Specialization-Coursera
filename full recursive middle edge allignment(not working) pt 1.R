# full recursive middle edge allignment

LinearSpaceAlignment <- function(v_raw, w_raw, match = 1, mismatch = 1, indel = 1) {
  v <- unlist(strsplit(v_raw, ""))
  w <- unlist(strsplit(w_raw, ""))
  if (length(v) > length(w)) {
    temp <- v
    v <- w
    w <- temp
  }
  
  alignment_edges <- list()
  
  Find_Middle_Edge <- function(v, w, match, mismatch, indel) {
    if (length(v) == 0 || length(w) == 0) {
      return(list(from = c(0, 0), to = c(0, 0), move = "→"))
    }
    n <- length(v)
    m <- length(w)
    mid <- floor(m / 2)
    score_fn <- function(a, b) ifelse(a == b, match, -mismatch)
    
    forward_dp <- function(v, w_sub) {
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
    
    reverse_dp <- function(v, w_sub) {
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
      return(rev(prev))
    }
    
    scoreL <- forward_dp(v, w[1:mid])
    scoreR <- reverse_dp(v, w[(mid + 1):m])
    total <- scoreL + scoreR
    i <- which.max(total) - 1
    
    diag <- if (i < n && mid < m) score_fn(v[i + 1], w[mid + 1]) else -Inf
    down <- -indel
    right <- -indel
    
    diag_score <- scoreL[i + 1] + diag + scoreR[i + 2]
    down_score <- if (i + 2 <= length(scoreL)) scoreL[i + 2] + down + scoreR[i + 2] else -Inf
    right_score <- scoreL[i + 1] + right + scoreR[i + 1]
    move <- which.max(c(diag_score, down_score, right_score))
    
    edge_moves <- list(
      c(i, mid, i + 1, mid + 1),
      c(i, mid, i + 1, mid),
      c(i, mid, i, mid + 1)
    )
    
    return(list(
      from = edge_moves[[move]][1:2],
      to   = edge_moves[[move]][3:4],
      move = c("↘", "↓", "→")[move]
    ))
  }
  
  build_alignment <- function(v, w, top, bottom, left, right) {
    if ((bottom - top) <= 0 || (right - left) <= 0) return()
    
    if ((bottom - top == 1) && (right - left == 1)) {
      move_type <- if (v[top + 1] == w[left + 1]) "↘" else "↘"
      alignment_edges[[length(alignment_edges) + 1]] <<- list(
        from = c(top, left),
        to   = c(top + 1, left + 1),
        move = move_type
      )
      return()
    }
    
    if (left == right) {
      for (i in top:(bottom - 1)) {
        alignment_edges[[length(alignment_edges) + 1]] <<- list(
          from = c(i, left),
          to   = c(i + 1, left),
          move = "↓"
        )
      }
      return()
    }
    
    if (top == bottom) {
      for (j in left:(right - 1)) {
        alignment_edges[[length(alignment_edges) + 1]] <<- list(
          from = c(top, j),
          to   = c(top, j + 1),
          move = "→"
        )
      }
      return()
    }
    
    mid <- floor((left + right) / 2)
    
    v_sub <- if (bottom - top > 0) v[top:(bottom - 1)] else character(0)
    w_sub <- if (right - left > 0) w[left:(right - 1)] else character(0)
    
    midEdge <- Find_Middle_Edge(v_sub, w_sub, match, mismatch, indel)
    from <- midEdge$from + c(top, left)
    to   <- midEdge$to + c(top, left)
    
    if (all(from == to)) return()
    
    build_alignment(v, w, top, from[1], left, from[2])
    alignment_edges[[length(alignment_edges) + 1]] <<- list(from = from, to = to, move = midEdge$move)
    build_alignment(v, w, to[1], bottom, to[2], right)
  }
  
  build_alignment(v, w, 0, length(v), 0, length(w))
  
  # Alignment formatter
  v_align <- ""
  w_align <- ""
  score <- 0
  for (edge in alignment_edges) {
    from <- edge$from
    move <- edge$move
    if (move == "↘") {
      vc <- v[from[1] + 1]
      wc <- w[from[2] + 1]
      v_align <- paste0(v_align, vc)
      w_align <- paste0(w_align, wc)
      score <- score + ifelse(vc == wc, match, -mismatch)
    } else if (move == "→") {
      v_align <- paste0(v_align, "-")
      w_align <- paste0(w_align, w[from[2] + 1])
      score <- score - indel
    } else if (move == "↓") {
      v_align <- paste0(v_align, v[from[1] + 1])
      w_align <- paste0(w_align, "-")
      score <- score - indel
    }
  }
  
  return(list(score = score, v_aligned = v_align, w_aligned = w_align, path = alignment_edges))
}


v <- "TGGCCTTAGTCAGGCTAGTACCTCATTGCGGATCCTGGGGAAACTTCATGGTGCATGCGGTGTATAATAGGTGCTCTCGCTCCACTATAATTACTGATCCCGTCCAAGAGGCCTAAGCCGCGTTACCTATAAAGGGCAACAGTCGCTACTAGTTAGCAGAAGGCGTCGGGTTCGACGGTCACACGGCTTGGCAACTCAGCTATGTGATCTTGACTAATGCATATGCAGCTGACTGAGGAAGAGTGGGGGTGGAATGCCTACCAACTTTAGTGGCCCCGTGTGTATGAGGATGAACTGGATGTAGAGATCGGTCAGAGCACGGCTTAGGAATACTATGGCCTTGATAACAAGGGTCTGGCGATAGAGCGAACGGAAGTGAACAGGTTGCTACTTCGTCGGACTGTTTTTTGCAATGCCATACAAAACCCTCTAAATTTAGGTTATGGAAATTGAAGATTATGTAACGTCGTCATACACGGCGTTCTTGATGCCATTTCTCGAGTGTGCTTCAAGTCCGGGGTTGATGGTTACGGTAAAAGCCTATGACCGAACCTTCAGGTCCAAACCCGGGGAATGCACTTGCCCAGCATCCTGAATAGGTCGACTTGTCCTCGACAAAATACGCTTGCAATGAGAATGGGACGTCATCACAGTTAATGATGCTTAAGCTCGCAGAGATATTACTGAGCTACGGCGGTAAACTCTCCGCAGTGTAAAACTCCAACACGTGAGGTAAGTTTGTATACAATCGACATCTCGCCTATGGACACTATGTGTCCCTCCGCGCGCTCCACTAGGGCGAAAATGAATCTGCTAGAGACCCTTGCTCATCAGACGAGTCGTTCGATATCGCATCCGACTCAGTGTTGAGTTCGTAAAGATAGGCTTGGTATACACCTACGACTATACCGTAAAACTTTCACCTCGGGCGGGTTGTGTTTGTGGAGGAAGTCCTTGTGGCGGGGAGCTCTGCCGCGAGTAGACGCGCGTTGCGGCCATCTAGCTCCTTCTACTGGCGTGGCCTCTTCAGAGCACCTCTTTGAGCAATAAAGGACCTTACTTGTAAGACGCTGCCGCGTGATGTTAGGCATAGCTCATGCAACCTCAGAGACCCTGCACAGCGAGTACAACTTGAATTACGCAGCACAGTTATCCGCAGGCCCGGGCGAGCAACAATTGAGACGACACTTATAGCCTAGCCTCGCCCTAATGACGTTGCTCGTGAGTCAAAGGACTTTACTCCACGGTTTCCTATGGACCAAAGTACGGCCTCGCGGCATAGATAGGCCATAGTGCGATAGCCCTTAGTGTTTCGGTTACGCTAAAACACACGTTGCATTCCCACGCATGGTAAAAAATTGGAGCAAATCTACCAACCGGCAAGTCTATCGTATAGGTAGTAGATCATAATCACCGCGAAAATTCGGTGCGAGTGCGGTGGAAGGGGTTAATACTTCCACTTATGTCCACCCCCTGCGTATTGCGATGCCAATTATCCTTGCTCTTCGGAGCAACAGGTCGCCATTTCCCCCGGTTTTGCAAAGCCGCCCAGGCCTGGGCAGCATGTTGGTATCCCCAGCGCTCAACGGTTCTCCTGGTGGCGCAGAGCACTGCTTGAAACAACAACGTGCGATTCCTGGGTCCGCGCGGTTCGTCCAGCTAGGTTTTAGAGTATGAGACATTATTGTTTGAAAGTAAACGGGTTCCGGACCCTTTCGGGCCTTTAGCCTGGAGGACTTTTAAAACGACCACTTGAACACTCCAAATCTTAGTTGAATCATTCTTGGCTGGTACAGGGAC"
w <- "TGGGGCTAGTACCGTGGAGAAACTTCAGCATGAATATGGCGTGCATGCCGTGCCGAGTGTATAATAGGTGCTCGGGCGCCGTCTAGCTAATTACTGTGGTATCACGTCCAAGAGGCCTAAACCACGTCTGTCATATCTATCTTGAAGGACCGGTTACCAACTGTCCGTCGGGTACTTCATCAGCGTTGTTATACCGCAGGGCCTCGGGTTCGACGGTCGGCAACTCAGCTATGTGCCAACTAACTGGCTCATATGCAGCTGACTGAGGATGCAGGGCGGCTTGGGGTGGAGCTGCGCGCTCCCAACCGACTCTAGAGGCCTCGTGTTCGGAGGGTGTGGCGTTCGGGGTTGAACTGAGAGGTCGAGATCGGTCGCTAGGAGATTAGTGATACAAGAAAAGAATGGCCTATGGTTTGATATGGACATCTCTGGGCGATACTTTAAAGCGAACGGACGTGTTATCTAGCCGCAGAGGTTCGTCGAAATAGACTGTTTTTCACAACGTATATGGCCATACAAAACCCTCTAAATTAGGTAACCTGGTTAAACCTGAGCTATGAAATTGATCCGGTACGTCGACGTCGCCTAGGCGTTGCTGCGATGCCATTTCTGACTAGCGGTGTGCTTCAAGTACAGGGTTGTACGGTAGAACTACAGGGATGGCACTTGCCCAGCAGGTCCACCTGGATAGGTCGACTTGTCCTCTGTGGTTCGACAAAATCCGCTTCCAATCCGGTGTGATAATGGGCCATCATCACAGTTGATGCTTAAGGCAGAGATATTGCTGAGCTCCCCCGGTACAATATCATCGGATCGTACCAACCCGTGAGGGAATACAATCGACACTCGATCGTGGATACGCCCCAAGTCCCTCCGCACGGGCCACTATATAGATGATCTGCTTGAGCCTTGCTCATCAGACGAGTGTTCGGTGTCGCTTCCCTCGGATGACGTTCATGATTGGTCGGCCAAGGTTCGTTGTTACCGCGTGGTAGACACCCTAGACCGTAAAACTAAGGCAGGTCCGTCCGTGGCGGCAGTTTGTGTTTGTAGGAATCATCCGCGCGGAGCTCTGCCGTTAAACTGCCGAGTAGACACGCGTTGCGGCCATCTAGGTGCTTCTACCGGATGGCCTGTGCTCCATCATTGCTCCGAGCAATACTAGAATTTCTTTGACTAATAAAGGCTGGACTGACCTTACTTGTAAGACGCTGCCGCGTGATGGTTACCGCAAACTCATGCAACCTTAGGCAGTTGCACAGCGAGTACAACTAGTTGAATTACGCAGGGCGGCAATGGATGTCTGGCTGAACTCGACGACAACGTGCCCCTCGCGATGACGTGGCTCGTGAGTCAAACGACTTTACTACACGGTTGTGAGTGGGCGGCAGATAGGCCATTGTCTAGACGCGATAGGTTACTTCTGTTGCAGTCTTCACGCATGGTCAAAGCAGTGCAAACAACCGTCTTTCCGATTAGAGGTCTATTGCAAATCACCGCGAAATTGTGCATCAGAGGGATCGGCGATGGCCTGGTTAATAATTCCACTTCTTCCCCAGTGTCGACCCCTGCGAAGTTTAGCATTGTATCCTTCAACAGTCGAATCTGCATTTCCCCAAAGCCATTGGCCACGCAGGGGGCAGCAGGTTGACATATCGATGTGGGGCCCGCGCTTTCTGGTGGCGGGTCCAGCCGTGAAACCCCGTGCAGGACCAATTGGTCCGCGAGCCAAACTACTCCAGCTAGGTTTTCGGACCGAGTAACATTTTGTAAACCCTTTCGGGTCGGGAGCTTGGAGGACTGTTAAAACGACCACTTGAACACTCACAAATCTCTCTCACTAAGTCTTGGCCTTTTTCGGGTACAGGGAC"


result <- LinearSpaceAlignment(v, w, 1, 1, 5)
# Open a connection to the file "results.txt"
file_conn <- file("results.txt")

# Create the content to write
output <- c(
  result$score,
  result$v_aligned,
 result$w_aligned)


# Write to the file
writeLines(output, file_conn)

# Close the connection
close(file_conn)
