### Full BWMatching

# BWT with sentinel appended
bwt_with_sentinel <- function(text) {
  text <- paste0(text, "$")
  n <- nchar(text)
  suffixes <- sapply(1:n, function(i) substr(text, i, n))
  sa <- order(suffixes)
  bwt_chars <- sapply(sa, function(pos) {
    if (pos == 1) substr(text, n, n) else substr(text, pos - 1, pos - 1)
  })
  paste(bwt_chars, collapse = "")
}

# Build FM-index components: First column, C array, and Occ prefix counts
build_bwt_index <- function(bwt_str) {
  last <- strsplit(bwt_str, "")[[1]]
  n <- length(last)
  first <- sort(last)
  
  # C[c] = first position (1-based) of character c in first column
  chars <- sort(unique(last))
  C <- setNames(sapply(chars, function(c) which(first == c)[1]), chars)
  
  # Occ[c][i] = # of occurrences of c in last[1:i]
  Occ <- lapply(chars, function(c) {
    # cumulative sum over last == c
    x <- as.integer(last == c)
    cs <- cumsum(x)
    # pad with 0 at 0-position for Occ(c, 0)
    c(0, cs)  # length n+1, Occ(c, i) accessed at index i
  })
  names(Occ) <- chars
  
  list(last = last, first = first, C = C, Occ = Occ, n = n, chars = chars)
}


BWMatching <- function(index, pattern) {
  last <- index$last
  C <- index$C
  Occ <- index$Occ
  n <- index$n
  
  # Work with character vector of pattern
  p <- strsplit(pattern, "")[[1]]
  
  top <- 1
  bottom <- n
  
  while (top <= bottom) {
    if (length(p) > 0) {
      symbol <- tail(p, 1)
      p <- head(p, -1)
      
      if (!(symbol %in% names(C))) {
        return(0)
      }
      
      # Occ access: Occ[[symbol]][i] where i is 0..n
      top_prime <- C[[symbol]] + Occ[[symbol]][top - 1 + 1] + 1
      bottom_prime <- C[[symbol]] + Occ[[symbol]][bottom + 1]
      
      if (top_prime > bottom_prime) {
        return(0)
      } else {
        top <- top_prime
        bottom <- bottom_prime
      }
    } else {
      return(bottom - top + 1)
    }
  }
  0
}

patterns <- "CCT CAC GAG CAG ATC"
patterns <- unlist(strsplit(patterns," "))
text <- "TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC"
bwt_str <- bwt_with_sentinel(text)
# bwt_str should be "annb$aa"
idx <- build_bwt_index(bwt_str)

result <- list()

for (i in seq_along(patterns)) {
  
  res <- BWMatching(idx, i)   
  result[[i]] <- res
}

###test

BWMatching <- function(bwt_str, patterns) {
  # Split BWT into vector of characters
  last <- strsplit(bwt_str, "")[[1]]
  n <- length(last)
  
  # First column = sorted BWT
  first <- sort(last)
  
  # Occurrence indices
  last_occ  <- ave(seq_along(last),  last,  FUN = seq_along)
  first_occ <- ave(seq_along(first), first, FUN = seq_along)
  
  # Last-to-First mapping
  LF <- match(paste(last, last_occ), paste(first, first_occ))
  
  # Function to count matches for one pattern
  count_pattern <- function(pattern) {
    pat <- strsplit(pattern, "")[[1]]
    top <- 1
    bottom <- n
    
    while (top <= bottom) {
      if (length(pat) > 0) {
        symbol <- tail(pat, 1)
        pat <- head(pat, -1)
        
        # positions from top..bottom in LastColumn
        positions <- top:bottom
        matches <- which(last[positions] == symbol)
        
        if (length(matches) == 0) {
          return(0)
        } else {
          topIndex <- positions[min(matches)]
          bottomIndex <- positions[max(matches)]
          top <- LF[topIndex]
          bottom <- LF[bottomIndex]
        }
      } else {
        return(bottom - top + 1)
      }
    }
    return(0)
  }
  
  # Apply to all patterns
  sapply(patterns, count_pattern)
}

# Example with your debug dataset
bwt_input <- "TTGACGCCCGCAGCTCTATCATAAGCCAAGCAATATATAGTCCTTCGATGGAGCCCTACAGCCCCACCGCTGTTGTGTTTTGGCGGAGATTGGGATACATTTACGGTAGAGGTGTCCGCCAAACGCGCCGCTCCTTCTCGGGTGTTGTACGTAACTCGACCATTTAAACCCAGGAGGTTCTTTGATCACGTGGCGTGGTGCAAGTACGGGATTGGGCTCCACCCTAGGCGATGCAACAAAGTTCCTAGTCTCCAAGACATGGGGATGATTTCCAGTATAAGCGGTGGCGTGCCCAGAGCTAGTAG$GGACATCACATCACCTAGGAACGTGTATCGTGTTCGAAGTCTCGGTAACCGTTAAATCGTGAAACTCCAGAGACGTAATATGTCTAAGGGTATCGGGGTTCGTATGCTGTTTAGTCCATTTTATCCCGTACTAGCAAGTAAGTGGGACACTAGACACCCACTTAGTCGTCCTTCTATACTAGTCTTCACTTTTTGACTATAACTTGGGTATCAGTTCATAGCCTGGGTTCCTTATAATGCACCAAAATTTAACCCTTCATCGTTGTTAACCGAGCATCGCCTCGGCATGTATCACCAGTGCCCCCGTACCTTGCTGCCCCACGCTTTAAGGCTGGTCTGAGTAATGGCGAAGAGCGCCCCAACCGTGATCTTCTAGAAAATCAGGCGACGCCGTGCTTAGGGGGACCTAATTTAGGGTCCCCTATCTGGTCAAGTAGGCCGAGGTCTCCGAGTCTATGGAGGGTGGTTAATCATGACAAGCGACCATTAGACTAAGGACTCCATGCACGACCCGTGGTCCTCTGAGCTCATACGCAGCCCTTCGAAACAAAACTAGGGAACGCCATAATCAGTGCAGTCATCGTGTCCTCCAGAGAAAGCCAACGTTTCGTGCAGAAAGAACCCCTTCGGTCGCTCCGGGAGCTTCCTCCATCGTTGATGGCCTTCATCATCAATGACATGTTCCCAACGAACTC"
input_file <- choose.files()
patterns <- readLines(input_file)
patterns <- unlist(strsplit(patterns," "))

input_file <- choose.files()
text <- readLines(input_file)

res <- BWMatching(bwt_input, patterns)
# Output: 2 1 1 0 1
cat(res)

