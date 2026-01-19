### Better BWMatching

# Build FM-index components: FirstOccurrence and Count
build_bwt_index <- function(bwt_str) {
  last <- strsplit(bwt_str, "")[[1]]
  n <- length(last)
  first <- sort(last)
  
  # FirstOccurrence[c] = first position (0-based) of character c in first column
  chars <- sort(unique(last))
  FirstOccurrence <- setNames(sapply(chars, function(c) {
    which(first == c)[1] - 1   # convert to 0-based
  }), chars)
  
  # Count[c][i] = # of occurrences of c in last[1:i], with Count(c,0)=0
  Count <- lapply(chars, function(c) {
    x <- as.integer(last == c)
    cs <- cumsum(x)
    c(0, cs)  # length n+1
  })
  names(Count) <- chars
  
  list(last = last, FirstOccurrence = FirstOccurrence, Count = Count, n = n)
}

BetterBWMatching <- function(index, pattern) {
  last <- index$last
  FirstOccurrence <- index$FirstOccurrence
  Count <- index$Count
  n <- index$n
  
  p <- strsplit(pattern, "")[[1]]
  top <- 0
  bottom <- n - 1
  
  while (top <= bottom) {
    if (length(p) > 0) {
      symbol <- tail(p, 1)
      p <- head(p, -1)
      
      if (!(symbol %in% names(FirstOccurrence))) {
        return(0)
      }
      
      top <- FirstOccurrence[[symbol]] + Count[[symbol]][top + 1]
      bottom <- FirstOccurrence[[symbol]] + Count[[symbol]][bottom + 1 + 1] - 1
      
      if (top > bottom) {
        return(0)
      }
    } else {
      return(bottom - top + 1)
    }
  }
  0
}

input_file <- choose.files()
full_text <- readLines(input_file)
full_text <- unlist(strsplit(full_text, "\n"))
patterns <- full_text[[2]]
patterns <- unlist(strsplit(patterns," "))

text <- full_text[[1]]

idx <- build_bwt_index(text)
cat(sapply(patterns, function(p) BetterBWMatching(idx, p))," ")
# Expected output: 1 2 1
