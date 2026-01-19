### Forward-Bacward HMM

soft_decode <- function(x, alphabet, states, transition_mat, emission_mat) {
  x <- strsplit(x, "")[[1]]
  n <- length(x); k <- length(states)
  
  # Forward
  forward <- matrix(0, nrow=n, ncol=k, dimnames=list(1:n, states))
  init <- rep(1/k, k)
  sym1 <- x[1]
  for (s in states) {
    forward[1,s] <- init[states==s] * emission_mat[s,sym1]
  }
  for (i in 2:n) {
    sym <- x[i]
    for (s in states) {
      forward[i,s] <- sum(forward[i-1,] * transition_mat[,s]) * emission_mat[s,sym]
    }
  }
  
  # Backward
  backward <- matrix(0, nrow=n, ncol=k, dimnames=list(1:n, states))
  backward[n,] <- 1
  for (i in (n-1):1) {
    sym <- x[i+1]
    for (s in states) {
      backward[i,s] <- sum(backward[i+1,] * transition_mat[s,] * emission_mat[,sym])
    }
  }
  
  # Sequence probability
  prob_x <- sum(forward[n,])
  
  # Posterior probabilities
  post <- matrix(0, nrow=n, ncol=k, dimnames=list(1:n, states))
  for (i in 1:n) {
    for (s in states) {
      post[i,s] <- (forward[i,s] * backward[i,s]) / prob_x
    }
  }
  
  post
}


# Pick file
input_file <- choose.files()

# Read and clean lines
lines <- readLines(input_file, warn = FALSE)
lines <- lines[nchar(trimws(lines)) > 0]  # drop empty lines

# Find separators "--------"
sep_idx <- which(grepl("^-{2,}$", trimws(lines)))

# 1) Emitted string (line before first separator)
x <- trimws(lines[1])

# 2) Alphabet (line after first separator)
alphabet <- unlist(strsplit(trimws(lines[sep_idx[1] + 1]), "\\s+"))

# 3) States (line after second separator)
states <- unlist(strsplit(trimws(lines[sep_idx[2] + 1]), "\\s+"))

# 4) Emission matrix (block after third separator)
mat_start <- sep_idx[3] + 1
mat_end <- sep_idx[4] - 1
mat_lines <- lines[mat_start:mat_end]

transition_mat <- as.matrix(read.table(text = mat_lines, header = TRUE, row.names = 1))

# 5) Emission matrix (block after 4th separator)
mat_start <- sep_idx[4] + 1
mat_lines <- lines[mat_start:length(lines)]

# Drop header row (alphabet labels)
mat_lines <- mat_lines[-1]

# Parse matrix: first column = row state labels, rest = numeric
mat_df <- read.table(text = mat_lines,
                     header = FALSE, stringsAsFactors = FALSE,
                     sep = "", fill = TRUE, strip.white = TRUE)

rownames_vec <- mat_df[[1]]
numeric_mat <- as.matrix(mat_df[, -1, drop = FALSE])
storage.mode(numeric_mat) <- "numeric"
rownames(numeric_mat) <- rownames_vec
colnames(numeric_mat) <- alphabet

emission_mat <- numeric_mat

post <- soft_decode(x, alphabet, states, transition_mat, emission_mat)
