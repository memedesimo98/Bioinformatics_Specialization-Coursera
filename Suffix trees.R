### Suffix trees

### Full suffix tree edges:

# Generate all suffixes of a string
suffixes <- function(text) {
  n <- nchar(text)
  sapply(1:n, function(i) substr(text, i, n))
}

# Build suffix tree edges by compressing common prefixes
suffix_tree_edges <- function(text) {
  sufs <- suffixes(text)
  
  # Recursive function to compress suffixes
  compress <- function(strings) {
    if (length(strings) == 0) return(character(0))
    if (length(strings) == 1) return(strings)
    
    # Group by first character
    groups <- split(strings, substring(strings, 1, 1))
    edges <- c()
    
    for (ch in names(groups)) {
      group <- groups[[ch]]
      
      # Find longest common prefix in this group
      lcp <- ch
      min_len <- min(nchar(group))
      k <- 1
      repeat {
        if (k >= min_len) break
        chars <- substring(group, k+1, k+1)
        if (length(unique(chars)) == 1) {
          lcp <- paste0(lcp, unique(chars))
          k <- k + 1
        } else break
      }
      
      # Add edge label
      edges <- c(edges, lcp)
      
      # Recurse on suffixes after the prefix
      rest <- substring(group, nchar(lcp)+1)
      rest <- rest[rest != ""]
      edges <- c(edges, compress(rest))
    }
    edges
  }
  
  compress(sufs)
}

# Example
text <- "TTGAATGACTCCTATAACGAACTTCGACATGGCA$"
edges <- suffix_tree_edges(text)
cat(edges)


### Most repeated suffix

# Longest Repeat Problem using suffix array + LCP
suffix_array <- function(text) {
  sufs <- suffixes(text)
  order(sufs)
}

lcp_array <- function(text, sa) {
  n <- length(sa)
  lcp <- integer(n)
  for (i in 2:n) {
    s1 <- substr(text, sa[i-1], nchar(text))
    s2 <- substr(text, sa[i], nchar(text))
    k <- 0
    while (k < min(nchar(s1), nchar(s2)) &&
           substr(s1, k+1, k+1) == substr(s2, k+1, k+1)) {
      k <- k + 1
    }
    lcp[i] <- k
  }
  lcp
}

longest_repeat <- function(text) {
  sa <- suffix_array(text)
  lcp <- lcp_array(text, sa)
  max_len <- max(lcp)
  if (max_len == 0) return("")
  idx <- which.max(lcp)
  start <- sa[idx]
  substr(text, start, start + max_len - 1)
}

# Solve Longest Repeat Problem
text2 <- "ATAGCTTACTGGAGAGCAGTGTTCATTTGGTATGTGCTAGACAGTCTTTTGTGGAGACAAACAAACTTCGGTACATCAAAATGGATTCGCAGCGTTACAAACGGAATAATGTCTGCTGATCCAGTAAGCAAAAGAGCTGCGTTTCACCGATTGCAGGGGGGATTCCGGCCACATCACGTGCGAGGGAATAGTAACGACCATACCAGTCCGACCCTATTTATTACTACCGCGAGCTCGGGTAAGAGTCTCTCTAACTGCTCCGAGTGTACTTAGCGGCAAGTAGATACACCCTCAGTCGTTAGTTTCTCTTACCTCCATCATGGCTACACAACTGTCCTAGTACTCGTGAACCCTCAAGCAAACAAACTTCAGCGCCCCATAGACCCTCGGGTAAGTCCTCGACAGTCGATGTCAACAATGCATCAATACTACAAAGCTTCCGGCCCTCGTATTGCCTCCCTTAGAGCAAGCCGCGCTTGTCCGTAGCTTTATTCGCTCCCCCACGCTGGGCCAAGGGTCCCTCCGCGTTACACGGCTGGCGACCCATTGGCTTTGCTCCAACGCGTCCGGTAACTCCAATATCGATAGGGAGTTTTTTGGTGTTCATCTTAGGGAATCAGCTATAGTCTTCCCCTGGTTACCTCAGATAGCTCTAGGATATCCCCTTAAATTATGCGCTGACTGCGTCACTCAATGACCGGGAAAGCTCTTTCAAGCACCGTATATGCTGATCATAATTGGTGGTTGGTAAGAGCAGAGATTTTGTGGAGACAAACAAACTTCGGTACATCAAAATGGATTCGCAGCGTTACAAACGGAATAATGTCTGCTGCGGACCGCTACAAAGGCGTCCATTGGACCACGACCGAGTCGTACTTTTAGAGAAGTTCGGCGTCTAGACTTACGAGAGACTCACTGTTTGTGGAGACAAACAAACTTCGGTACATCAAAATGGATTCGCAGCGTTACAAACGGAATAATGTCTGCTGACCCTCGCGCTGCTCCCCGACTTGTGCCCTACCAAATCCCCTCCCCAACTGTCAAGCAAGCGCATCCCTGTTTCTTGAGTGGACACAACAACGTAGTTTATAAGAGCTTGATCCCTATTCTGACTTGGAATGGGGGCCCACAATTTAGTATCCGGGCTTTAGAGTCTAGGTTTACTTCCTTCAAGCCTTACCCTGACCGTTACCCGTTGTGTCTGGTCCGCGGG"
repeated <- longest_repeat(text2)
print("Longest Repeat:")
print(repeated)

### most repeated suffix over 2 texts

# Generate all suffixes
suffixes <- function(text) {
  n <- nchar(text)
  sapply(1:n, function(i) substr(text, i, n))
}

# Build suffix array
suffix_array <- function(text) {
  sufs <- suffixes(text)
  order(sufs)
}

# Compute LCP array
lcp_array <- function(text, sa) {
  n <- length(sa)
  lcp <- integer(n)
  for (i in 2:n) {
    s1 <- substr(text, sa[i-1], nchar(text))
    s2 <- substr(text, sa[i], nchar(text))
    k <- 0
    while (k < min(nchar(s1), nchar(s2)) &&
           substr(s1, k+1, k+1) == substr(s2, k+1, k+1)) {
      k <- k + 1
    }
    lcp[i] <- k
  }
  lcp
}

# Longest Shared Substring
longest_shared_substring <- function(text1, text2) {
  # Concatenate with separator
  combined <- paste0(text1, "#", text2, "$")
  sa <- suffix_array(combined)
  lcp <- lcp_array(combined, sa)
  
  max_len <- 0
  result <- ""
  n1 <- nchar(text1)
  
  for (i in 2:length(sa)) {
    # Check if suffixes come from different strings
    from1 <- sa[i-1] <= n1
    from2 <- sa[i] <= n1
    if (from1 != from2) {
      if (lcp[i] > max_len) {
        max_len <- lcp[i]
        result <- substr(combined, sa[i], sa[i] + max_len - 1)
      }
    }
  }
  result
}

# -------------------------------
# Example usage
# -------------------------------
text <- "AAAATAAACAAAGAATTAATCAATGAACTAACCAACGAAGTAAGCAAGGATATACATAGATTTATTCATTGATCTATCCATCGATGTATGCATGGACACAGACTTACTCACTGACCTACCCACCGACGTACGCACGGAGAGTTAGTCAGTGAGCTAGCCAGCGAGGTAGGCAGGGTTTTCTTTGTTCCTTCGTTGCTTGGTCTCTGTCCCTCCGTCGCTCGGTGTGCCTGCGTGGCTGGGCCCCGCCGGCGCGGGGAAAAAATGCGGATTCGTGGAGCGGGCGTCCACTAGAACGTAGATGGCGGCACAGGCGAGGTTTCGGCCCGCGCCGTCGCGGTAACATGCTCCTCTCCTTTCATTGGTCTTAATAGGCTGTATAATTGAGTAGACTGTACTCCAATTACTAATGGAGTAGCATGGCAAGCGGTAACAGACACTAGACCTATGAGCCGTAACTCGCACACAACAGTAGCCCGGGCTCTCCGGGATCTACTCGGCTCCAATCGCCTGCTATCTCTAGCCGCTGTACGACAGGTCTCGACACTATAGATAAGTACGCTTGCTGCCGCCCAGGTTCACTTAATCAACTGGAGTAGGAAATATGTCGAAACGTAATGCCATTCAATGTACGTTGCACTCCCAATGTTCGAGGACCTGGG
AAAATAAACAAAGAATTAATCAATGAACTAACCAACGAAGTAAGCAAGGATATACATAGATTTATTCATTGATCTATCCATCGATGTATGCATGGACACAGACTTACTCACTGACCTACCCACCGACGTACGCACGGAGAGTTAGTCAGTGAGCTAGCCAGCGAGGTAGGCAGGGTTTTCTTTGTTCCTTCGTTGCTTGGTCTCTGTCCCTCCGTCGCTCGGTGTGCCTGCGTGGCTGGGCCCCGCCGGCGCGGGGAAA"
text1 <- unlist(strsplit(text,"\n")) [1]
text2 <- unlist(strsplit(text,"\n")) [2]

shared <- longest_shared_substring(text1, text2)
print(shared)

### shortest non shared

# Function to find the shortest non-shared substring
shortest_non_shared <- function(text1, text2) {
  n1 <- nchar(text1)
  
  # Check substrings by increasing length
  for (len in 1:n1) {
    for (i in 1:(n1 - len + 1)) {
      sub <- substr(text1, i, i + len - 1)
      if (!grepl(sub, text2, fixed = TRUE)) {
        return(sub)
      }
    }
  }
  return("")  # if all substrings are shared
}

# Example usage
text1 <- "CCAAGCTGCTAGAGG"
text2 <- "CATGCTGGGCTGGCT"

result <- shortest_non_shared(text1, text2)
print(result)
