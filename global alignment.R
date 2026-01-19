# LCS_backtracer_scoring

LCS_backtracer_scoring <- function(v, w, match, mismatch, indel) {
  v <- unlist(strsplit(v, "")) 
  w <- unlist(strsplit(w, "")) 
  V <- length(v) 
  W <- length(w) 
  s <- matrix(0, nrow = V+1, ncol = W+1) 
  Backtrack <- matrix("", nrow = V+1, ncol = W+1) 
  # Initialize score matrix 
  s[ ,1] <- 0:-V * indel 
  s[1, ] <- 0:-W * indel 
  # Initialize Backtrack for boundaries 
  Backtrack[2:(V+1), 1] <- "d" 
  # First column: deletions from v 
  Backtrack[1, 2:(W+1)] <- "r" 
  # First row: insertions from w 
  for (i in 2:(V+1)) {
    for (j in 2:(W+1)) {
      # Compute all three scores
      diag_score <- s[i - 1, j - 1] + ifelse(v[i - 1] == w[j - 1], match, -mismatch)
      up_score   <- s[i - 1, j] - indel
      left_score <- s[i, j - 1] - indel
      
      max_score <- max(diag_score, up_score, left_score)
      s[i, j] <- max_score
      
      # Prefer diagonal > up > left (you can adjust priority here)
      if (diag_score == max_score) {
        Backtrack[i, j] <- "dr"
      } else if (up_score == max_score) {
        Backtrack[i, j] <- "d"
      } else {
        Backtrack[i, j] <- "r"
      }
      
    } 
  } 
  alignment <- Output_Alignment(Backtrack,v,w,V+1,W+1)
  return(list(Backtrack = Backtrack, matrix = s,End_score = s[V+1,W+1], Alignment = alignment) )
  }

v <- "CGGCTCTGAATCTGTTACTGGTGACTACGTTTAGAGTTATGTGACTGGGCTTTGCTAGCTCGGATCTTCTCAGGGTTGAAACCGCTTTGAAACCTTGTGCCGAGTTGTCCCTCGGGGTCTACTACCCAGATTCGCGGGGGGACGAGTACGCCTTCTTGAGATTTATGTGCCTACTCCTACAAGGCACGCTTAGGTATGGACACCCCATCTGGGGTCTGCTGCATATGCAACCTGTCCTTCTTCGGGAACGTGTAAAGCGAAACTTATGTCCCATCCTCTTGTATGGAGAGCATGCGGATCTGCAGATATCAGATCCCCCACATCGTGCAACAGAATCGACCTTCTTTAGCGAAAAATATCTTTTGAGACGTAAGCTGAAGTCTCGGGAACTTGTCTACGCTTTAATTTCGACACACTCTTACCAGTGCCTTGGCAGCATCGGGATTCGTCGCAGCAGGGACAGCTTGCCAAATGAGAAACTCGAGACTAGTCTTCATACATATGTCCGGGTATTGTATGCCCAGATAGCAGAAGCCCAAAGATCGAGCTCTAGCGAGGAGCATAGCGTAAGCTCACCAATCCAAGACTGTAATGTTCCATTGCGTTACCGACATAACTAGAATTGGACGAACAGATTGGTTCGCGGCAGCATTGTCAACGGTGCATTTCCACCCCGTACAACTGACGGACAAATTACGCGTCAACGTAATGGGGGCCAGCCAAAAGATTGGGGCGTTAGTACGCAGAGTTTCGCCTCCAGGGGAGGCGTCGCTGGCGATCGTCGAGTATTAGGGCAGAGAGGCAAGTCCCCTACGCGACGCCGCCGGTCGAGTCTTCGAATTGGGAAAACACTGACGGGCAGTGTCAAGCAATTCTTGAACAGGCACCATCGCGTCCAGG"
w <- "AGTCGGGCCATCTGGGCTACGGTGACAAACAAGACGTTATACCCCTAGAGTTAAGATGTAATGGTTCAGGGTCTATAAGGGTTGTCACCCTGTTCTCTTGAGCTAGTCCCCAGTACAAAGTTGTCCCTCTACTGCACGGGGGGACGAGTACGCCTTCTTGAAATGCTTTATGTGCCTACTCCTGCGATGCAACCTGAAGCACCCATGCGCTTAGGAATGGACACCCCATCTGGGGTCTGCTGCATATGAGAGAGCTGCCCTAGGCGGGAACGTGTTAAGCGAAAATTCATGTCCCAGTGGCCCTGCCTTGTAAGCATGCGGATCTGCAGAAAACTAGATTAGATCTCCCACATCGTGCAGAACCGACAGTAGCTCTAACACTTGACCTTCTCTAGCGAAATCAACCCTGTCTTTTAAGACCCCGTTCAAGCTGAAGTCTCGGGAACAGGTCTACGCTTTAATTCCGGCACCACAGCTGGATTTGGCTGTCCAGGGCCTGCTTTCGTCGCAGCAGGTGAGAGTTAGTTTGCTAAAAGAGTAACGACTAACTCCGGGTACTTTCGTTTTGTACCAGATATCGTCCAATGATCAGCCGTCGAGATCTAGCGATTACCAGAGGGGCTTGGGTAGTGAAGTAGTGAGCTAGCCAAACCAATCCAAGACTTATGAGGCGGGTCCGCTCCATACTTCTCGCTAGAATTGGATGGAAAAGTGCGAGTTGCGATTGGTTCGCGGCAGCATAGTCAACGGTGCAATTCCACCCCGTAAGACCAGTTAGTGGCACGGACAGATTAACTTTTATGCCTGATAAACTCAACCTTAATGGGGGCCAAAGATTGGGGCTGACGTTAGTACGTTCGGAGTTTTGCCTAGGGGAGGCGTGAGTATTAGCGCAGGAGGCAACTCCCCTACGCGATGCCGCCGGTCGAGTCTTCGAATTGGGAAAACACTGACGGGCAGTCCGTCAAGCTACTCTTTTTAACAGGCCCCCGCGTCCAGG"
match <- 1
mismatch <- 1
indel <- 5
result <- LCS_backtracer_scoring(v,w,match,mismatch,indel)

Output_Alignment <- function(Backtrack, v, w, i, j) {
  if (i == 1 && j == 1) {
    return(list(v_align = "", w_align = ""))
  }
  
  if (i > 1 && Backtrack[i, j] == "d") {
    res <- Output_Alignment(Backtrack, v, w, i - 1, j)
    return(list(v_align = paste(res$v_align,v[i - 1],sep=""),
                w_align = paste( res$w_align,"-",sep="")))
  }
  
  if (j > 1 && Backtrack[i, j] == "r") {
    res <- Output_Alignment(Backtrack, v, w, i, j - 1)
    return(list(v_align = paste( res$v_align,"-",sep=""),
                w_align = paste( res$w_align,w[j - 1],sep="")))
  }
  
  if (i > 1 && j > 1 && Backtrack[i, j] == "dr") {
    res <- Output_Alignment(Backtrack, v, w, i - 1, j - 1)
    return(list(v_align = paste( res$v_align,v[i - 1],sep=""),
                w_align = paste( res$w_align,w[j - 1],sep="")))
  }
  
  return(list(v_align = "", w_align = ""))
}
result$End_score
result$Alignment

