# overlap alignment

# Backtrack function: reconstructs the optimal alignment path
Output_Alignment <- function(Backtrack, v, w, i, j) {
  if (i == 1 && j == 1) {
    return(list(v_align = "", w_align = ""))
  }
  
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

# Main Overlap Alignment function
Overlap_Alignment <- function(v, w, match = 1, mismatch = 1, indel = 2) {
  v <- unlist(strsplit(v, "")) 
  w <- unlist(strsplit(w, "")) 
  V <- length(v) 
  W <- length(w) 
  
  # Initialize score and backtrack matrices
  s <- matrix(0, nrow = V+1, ncol = W+1) 
  Backtrack <- matrix("", nrow = V+1, ncol = W+1) 
  
  # Penalize unmatched prefixes of both v and w (global-style init)
  s[, 1] <- 0:-V * indel
  s[1, ] <- 0:-W * indel
  
  Backtrack[2:(V+1), 1] <- "d"
  Backtrack[1, 2:(W+1)] <- "r"
  
  # Fill in dynamic programming matrix
  for (i in 2:(V+1)) {
    for (j in 2:(W+1)) {
      diag_score <- s[i - 1, j - 1] + ifelse(v[i - 1] == w[j - 1], match, -mismatch)
      up_score   <- s[i - 1, j] - indel
      left_score <- s[i, j - 1] - indel
      
      max_score <- max(diag_score, up_score, left_score)
      s[i, j] <- max_score
      
      if (left_score == max_score) {
        Backtrack[i, j] <- "r"  # insertion (gap in v)
      } else if (up_score == max_score) {
        Backtrack[i, j] <- "d"  # deletion (gap in w)
      } else {
        Backtrack[i, j] <- "dr" # match/mismatch
      }
    } 
  }
  
  # 🧠 Find the max score for any suffix of v vs any prefix of w
  best_score <- -Inf
  best_i <- 1
  best_j <- 1
  for (i in 2:(V+1)) {
    for (j in 2:(W+1)) {
      if (s[i, j] > best_score) {
        best_score <- s[i, j]
        best_i <- i
        best_j <- j
      }
    }
  }
  
  alignment <- Output_Alignment(Backtrack, v, w, best_i, best_j)
  
  return(list(
    Score = best_score,
    Alignment = alignment,
    AlignmentEndsAt = list(i = best_i, j = best_j),
    Matrix = s,
    Backtrack = Backtrack
  ))
}


v <- "CGGGATTGGAAGAAGGGTACCGCCCAATACCAGACAACCCACAAATTGAGCGCACAACTTCTTAACCGAGATGCAAATTGAGAAACGGACCAACAACCTGTACGTTGTTCCTGTTTCTATGAATAAGTCCACCTTGCAGATTTATTACAATAGCCTCGAATGCGCACCCAGATGTTTGGCCAAATCCGGAGTGGCGGTAGTCGTTGGGGACGTGCTTGTACACCCCACACAACGGATTTCATCTGGGGTTTCATCGATATGGGAGCGCTTAGGTAGGAACATCTCAAGGAGTGTCCGTCATAGATGTAATATGCCTAATTGTCACCATGCTCTTCGTACCCCCAGCCACGGAATGCTGGATTCTCAGAGGTGTGTCACTACTAGTAGCCCCCAATGATACGATGTTTTGCCAAAAAGGGGAAACAAGCCATCATGAGGTCAGATCGGTCGTTAGGTCGGAAAAGTGAGGTCTTTAGAGATACCCGCACGTAGTTTAGCCTCCATTCTCTGGGAGGGTCATGCTACTAAATATCGTCGTGGTCTGAATGCGGGCGACGCGTTATTAATAAATTTCATAATCAACGACAGCGGCATAACACTCACATATACCGGTGGACTTCGGTCCTGATCCACGTAGCAGTCTATCCCCCCATAATTCTAGCAGCGAGTCCACAAAGGGGATAATTGTTTAACGTCTCGATAAACGACCGCTTATGCGTTTCAGTCGGTATTACCCATAACACTAAGTGGCGGAACTCTTGACACGTCTATGAGCTTAGGCGTGCATGTTCGAGATCAGTCAACGTAACGCCGTATGGAGGCGCACGTATAATCATAATGGCGGGCGAAGCACCCCTTGCAGCATAGATATGGAGGTAAACTCGCCCCCATAGCCTACGA"
w <- "CGCTAAAAAAAACTTGGAAGACGGGTACCGCCCAATACCAGATTCTAGCCATGAATATTTAAGCACAAAGCCAATAACGGAGATAGAGTACAAATTGAGATATAACGGACCAACAACCTGTGACCTACGTTGTTAGCCTGTATCTATGAATACTTGCAGATTTATTAGGGTTAAATAGATCGAATGCCCTTTGGCCAAATCCTATACCCGGATTGGCGGTAGTCGATTGGTTGGTACACCCCACATATCGGTGGGGTTCGACAATGGGAGCGCTCAGGTAGGAACCCATAGAAGCCTCAAGGACTGCATAGATGTAATATGCCGACTGAGAATTGTCTCCATGCTCTTCGTACCTGACCAGTGACGGAATCGAGTGCTTGAACTAGTCTCAGTTGACAAAGAGGAGGAGGTATTCAACTACTAGTGACCCAGAGGCCCCAAGAGCGTTGCATGATACGATGTTTTGGCAAGGGAAACAAGCCATAAGAGCGCGAATAAGTCCGCACAAGAACGCAGTCCTTAGCCTCCATTCTATGGGAGGACTAAACACTGGGCTTGAATGGACGTTGTCGTATGAATTGCATAGAAGTGACAGCGAGCTAAATATAACACATACATATAGCCGTGGTACGATCTTCTCAGCCTGATCCACGTAGCAGTCATGCTGAATTCTAGCAGCGAGAGAAGGGATTGTTCTAACGTCTCGATGAACGACCGCCGCCTATGTACGCGTTTCAGTCGGTTGCGGATAATTTTAATGTACTAAGTAGCGGAACTCCGTTATGAGCTCCTGACGGTGCATGTTCGAGATTAGTCAGCGTAGTACGGAGGCGCTATAATCATAATGGCGAGCGAAGCACCCCTCTATTGAGTGTAACCTCGGTTGATACCCGCATAGCCTCCGA"

# Suppose 'result' is the output of Overlap_Alignment(...)
result <- Overlap_Alignment(v, w)

# Write alignment to text file
output_text <- paste(
  result$Score,
  result$Alignment$v_align,
  result$Alignment$w_align,
  sep = "\n"
)

writeLines(output_text, "alignment_output.txt")

Overlap_Alignment(v, w)
