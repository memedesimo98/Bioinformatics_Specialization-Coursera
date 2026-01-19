### suffix array and BWT


suffixes_with_positions <- function(text) {
  n <- nchar(text)
  # collect suffixes with their starting positions
  sfx <- sapply(1:n, function(i) substr(text, i, n))
  # get the order of suffixes alphabetically
  ord <- order(sfx)
  # return both sorted suffixes and their positions
  data.frame(
    position = ord - 1,        # starting positions in the original text
    suffix   = sfx[ord]    # suffixes in sorted order
  )
}


text <- "banana$"
res <- suffixes_with_positions(text)
cat(res$position, " ")

bwt <- function(text) {
  n <- nchar(text)
  
  # naive suffix array (O(n log n))
  suffixes <- sapply(1:n, function(i) substr(text, i, n))
  sa <- order(suffixes)
  
  # build BWT
  bwt_chars <- sapply(sa, function(pos) {
    if (pos == 1) substr(text, n, n) else substr(text, pos - 1, pos - 1)
  })
  
  paste(bwt_chars, collapse = "")
}

bwt("CGTTTGCTAT$")
# Output: "annb$aa"

inverse_bwt_df <- function(bwt_str) {
  n <- nchar(bwt_str)
  bwt_vec <- strsplit(bwt_str, "")[[1]]
  
  # First column = sorted BWT
  first_col <- sort(bwt_vec)
  
  # Occurrence indices
  bwt_occ   <- ave(seq_along(bwt_vec), bwt_vec, FUN = seq_along)
  first_occ <- ave(seq_along(first_col), first_col, FUN = seq_along)
  
  # Build data frame
  df <- data.frame(
    bwt       = bwt_vec,
    bwt_occ   = bwt_occ,
    first     = first_col,
    first_occ = first_occ,
    stringsAsFactors = FALSE
  )
  
  # LF-mapping
  df$lf <- match(paste(df$bwt, df$bwt_occ),
                 paste(df$first, df$first_occ))
  
  # Reconstruct original text
  res <- character(n)
  idx <- which(df$bwt == "$")  # start at sentinel
  for (i in n:1) {             # walk backwards to fill string forwards
    res[i] <- df$bwt[idx]
    idx <- df$lf[idx]
  }
  
  paste(res, collapse = "")
}

# Example
inverse_bwt_df("AT$AAACTTCG")
# "banana$"


### partial suffix array

partial_suffix_array <- function(text, K) {
  n <- nchar(text)
  # build full suffix array
  suffixes <- sapply(0:(n-1), function(i) substr(text, i+1, n))
  sa <- order(suffixes) - 1  # 0-based positions
  
  # filter for divisibility by K
  res <- list()
  for (i in seq_along(sa)) {
    if (sa[i] %% K == 0) {
      res[[length(res)+1]] <- paste(i-1, sa[i])  # i-1 for 0-based index
    }
  }
  return(res)
}

# Example
input_file <- choose.files()
text <- readLines(input_file)
K <- as.integer(text[[2]])
text <- text[[1]]
if (!grepl("\\$", text)) {
  text <- paste0(text, "$")
}
res <- unlist(partial_suffix_array(text, K))
writeLines(res, "result.txt")

