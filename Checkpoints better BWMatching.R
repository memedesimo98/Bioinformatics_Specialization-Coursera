### Checkpoints better BWMatching

# BWT with sentinel
bwt_with_sentinel <- function(text) {
  if (!grepl("\\$", text)) {
    text <- paste0(text, "$")
  }
  n <- nchar(text)
  suffixes <- sapply(1:n, function(i) substr(text, i, n))
  sa <- order(suffixes)
  bwt_chars <- sapply(sa, function(pos) {
    if (pos == 1) substr(text, n, n) else substr(text, pos - 1, pos - 1)
  })
  paste(bwt_chars, collapse = "")
}


# Count(symbol, i) using checkpoints (0 <= i <= n)
count_symbol <- function(index, symbol, i) {
  if (i == 0) return(0L)
  last <- index$last
  checkpoints <- index$checkpoints[[symbol]]
  all_idx <- as.integer(names(checkpoints))
  nearest <- max(all_idx[all_idx <= i])
  count <- checkpoints[as.character(nearest)]
  if (nearest < i) {
    count <- count + sum(last[(nearest + 1):i] == symbol)
  }
  as.integer(count)
}

# BetterBWMatching with checkpoints: returns first-column interval [top..bottom] (0-based)
BetterBWMatchingCheckpoints <- function(index, pattern) {
  p <- strsplit(pattern, "")[[1]]
  top <- 0L
  bottom <- index$n - 1L
  
  while (top <= bottom) {
    if (length(p) > 0) {
      symbol <- tail(p, 1)
      p <- head(p, -1)
      if (!(symbol %in% names(index$FirstOccurrence))) return(integer(0))
      
      # Correct 0-based update
      top    <- index$FirstOccurrence[[symbol]] + count_symbol(index, symbol, top)
      bottom <- index$FirstOccurrence[[symbol]] + count_symbol(index, symbol, bottom + 1L) - 1L
      
      if (is.na(top) || is.na(bottom) || top > bottom) return(integer(0))
    } else {
      return(top:bottom)
    }
  }
  integer(0)
}

# Solve multiple patterns: map first-column interval to last-column rows, then to positions
solve_multiple_patterns <- function(text, patterns, C = 100) {
  index <- build_index_with_checkpoints(text, C)
  index <- compute_distances(index)  # one-time O(n) preprocessing
  results <- lapply(patterns, function(pat) {
    interval <- BetterBWMatchingCheckpoints(index, pat)
    if (length(interval) == 0) return(paste0(pat, ":"))
    rows_last <- index$invLF[interval + 1L]
    pos <- (index$n - 1L) - index$dist_to_dollar[rows_last]
    pos <- sort(pos)
    paste0(pat, ": ", paste(pos, collapse = " "))
  })
  cat(paste(results, collapse = "\n"), "\n")
  return(results)
}



#### print test

compute_distances <- function(index) {
  n <- index$n
  dist <- rep(NA_integer_, n)
  dist[index$dollar_row_last] <- 0L
  for (r in seq_len(n)) {
    if (is.na(dist[r])) {
      path <- integer(0)
      cur <- r
      cat("Starting row:", r, "\n")
      while (is.na(dist[cur])) {
        cat("  Visiting row:", cur, "next_last:", index$next_last[cur], "\n")
        path <- c(path, cur)
        cur <- index$next_last[cur]
        # safety check
        if (length(path) > n) {
          stop("Loop detected in compute_distances at row ", r)
        }
      }
      base <- dist[cur]
      for (k in seq_along(path)) {
        dist[path[k]] <- base + k
        cat("    Assigning dist[", path[k], "] =", base + k, "\n")
      }
    }
  }
  index$dist_to_dollar <- dist
  index
}

build_index_with_checkpoints <- function(text, C = 100) {
  bwt <- bwt_with_sentinel(text)
  last <- strsplit(bwt, "")[[1]]
  n <- length(last)
  
  # First column (only to build FirstOccurrence and LF)
  first <- sort(last)
  
  # FirstOccurrence (0-based)
  chars <- sort(unique(last))
  FirstOccurrence <- setNames(sapply(chars, function(c) {
    which(first == c)[1] - 1
  }), chars)
  
  # Checkpoints for Count (ensure index 0 exists)
  checkpoints <- lapply(chars, function(c) {
    x <- as.integer(last == c)
    cs <- c(0, cumsum(x))   # indices 0..n
    idx <- unique(c(0, seq(0, n, by = C)))
    vals <- cs[idx + 1]
    names(vals) <- idx
    vals
  })
  names(checkpoints) <- chars
  
  # LF and inverse LF (1-based)
  last_rank  <- ave(seq_along(last),  last,  FUN = seq_along)
  first_rank <- ave(seq_along(first), first, FUN = seq_along)
  LF    <- match(paste(last,  last_rank),  paste(first, first_rank))  # last -> first
  invLF <- match(paste(first, first_rank), paste(last,  last_rank))   # first -> last
  
  # Next step in last-column space: follow LF directly
  next_last <- LF
  
  # Distance to sentinel row in last-column space (0-based positions)
  dollar_row <- which(last == "$")[1]
  dist_to_dollar <- rep(NA_integer_, n)
  dist_to_dollar[dollar_row] <- 0L
  
  # Fill distances in O(n) via memoized traversal
  for (r in seq_len(n)) {
    if (is.na(dist_to_dollar[r])) {
      path <- integer(0)
      cur <- r
      cat("Starting row:", r, "\n")
      while (is.na(dist_to_dollar[cur])) {
        cat("  Visiting row:", cur, "next_last:", next_last[cur], "\n")
        path <- c(path, cur)
        cur <- next_last[cur]
        if (length(path) > n) {
          stop("Loop detected in build_index_with_checkpoints at row ", r)
        }
      }
      base <- dist_to_dollar[cur]
      for (k in seq_along(path)) {
        dist_to_dollar[path[k]] <- base + k
        cat("    Assigning dist[", path[k], "] =", base + k, "\n")
      }
    }
  }
  
  list(
    bwt = bwt, last = last, FirstOccurrence = FirstOccurrence,
    checkpoints = checkpoints, C = C, n = n,
    LF = LF, invLF = invLF,
    next_last = next_last,
    dollar_row_last = dollar_row,
    dist_to_dollar = dist_to_dollar
  )
}

# Example
text <- "ATATATATAT"
patterns <- "GT AGCT TAA AAT AATAT"
patterns <- unlist(strsplit(patterns," "))
system.time( res <- solve_multiple_patterns(text, patterns,C = 100))
writeLines(unlist(res), "result.txt")


# Output:
# ATCG: 1 11
# GGGT: 4 15

