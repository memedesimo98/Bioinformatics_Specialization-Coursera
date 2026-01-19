### HMM decoding problem (Viterbi)

viterbi_decode <- function(x, alphabet, states, transition_mat, emission_mat) {
  x <- strsplit(x, "")[[1]]
  
  # Checks
  stopifnot(all(x %in% alphabet))
  stopifnot(all(rownames(transition_mat) == states),
            all(colnames(transition_mat) == states))
  stopifnot(all(rownames(emission_mat) == states),
            all(colnames(emission_mat) == alphabet))
  
  n <- length(x)
  k <- length(states)
  state_idx <- setNames(seq_len(k), states)
  
  logT <- log(transition_mat)
  logE <- log(emission_mat)
  logInit <- rep(log(1 / k), k)
  
  dp <- matrix(-Inf, nrow = n, ncol = k)
  colnames(dp) <- states
  bp <- matrix(NA_integer_, nrow = n, ncol = k)
  colnames(bp) <- states
  
  # Initialization
  sym1 <- x[1]
  for (s in states) {
    si <- state_idx[[s]]
    dp[1, si] <- logInit[si] + logE[s, sym1]
  }
  
  # Recurrence
  for (i in 2:n) {
    sym <- x[i]
    for (s in states) {
      si <- state_idx[[s]]
      candidates <- dp[i - 1, ] + logT[, s]
      best_prev_idx <- which.max(candidates)
      dp[i, si] <- candidates[best_prev_idx] + logE[s, sym]
      bp[i, si] <- best_prev_idx
    }
  }
  
  # Termination
  last_best_idx <- which.max(dp[n, ])
  best_log_prob <- dp[n, last_best_idx]
  
  best_path_idx <- integer(n)
  best_path_idx[n] <- last_best_idx
  for (i in (n - 1):1) {
    best_path_idx[i] <- bp[i + 1, best_path_idx[i + 1]]
  }
  
  best_states <- states[best_path_idx]
  best_path <- paste(best_states, collapse = "")
  
  # Return both path and probability
  list(
    path = best_path,
    log_prob = best_log_prob,
    prob = exp(best_log_prob)
  )
}


# Sample input components
x_str <- "xyxzzxyxyy"
alphabet <- c("x","y","z")
states <- c("A","B")
transition_mat <- matrix(
  c(0.303, 0.697,
    0.831, 0.169),
  nrow = 2, byrow = TRUE,
  dimnames = list(states, states)
)
emission_mat <- matrix(
  c(0.533, 0.065, 0.402,
    0.342, 0.334, 0.324),
  nrow = 2, byrow = TRUE,
  dimnames = list(states, alphabet)
)


# Run
decoded <- viterbi_decode(x_str, alphabet, states, transition_mat, emission_mat)
decoded$path
decoded$log_prob
# Expected: "AAABBAAAAA"

# Pick file
input_file <- choose.files()

# Read and clean lines
lines <- readLines(input_file, warn = FALSE)
lines <- lines[nchar(trimws(lines)) > 0]  # drop empty lines

# Find separators "--------"
sep_idx <- which(grepl("^-{2,}$", trimws(lines)))

# 1) Emitted string (line before first separator)
string <- trimws(lines[1])

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

decoded <- viterbi_decode(lines, alphabet, states, transition_mat, emission_mat)
#### test

states <- c("I","N")

transition_mat <- matrix(
  c(0.999,  0.001,   # I->I, I->N
    0.0001, 0.9999), # N->I, N->N
  nrow = 2, byrow = TRUE,
  dimnames = list(states, states)
)

alphabet <- c("A","C","G","T")

emission_mat <- matrix(
  c(0.10, 0.40, 0.40, 0.10,   # I emissions
    0.30, 0.20, 0.20, 0.30),  # N emissions
  nrow = 2, byrow = TRUE,
  dimnames = list(states, alphabet)
)

# Pick file
input_file <- choose.files()

# Read lines from FASTA file
lines <- readLines(input_file)

# Remove header lines (those starting with ">")
seq_lines <- lines[!grepl("^>", lines)]

# Collapse into one long string
seq_raw <- paste(seq_lines, collapse = "")

# Normalize case
seq_upper <- toupper(seq_raw)

# Keep only A, C, G, T
seq_clean <- gsub("[^ACGT]", "", seq_upper)



decoded <- viterbi_decode(seq_clean, alphabet, states, transition_mat, emission_mat)
decoded

count_islands <- function(path_str, island_state = "I") {
  s <- strsplit(path_str, "")[[1]]
  rle_states <- rle(s)
  sum(rle_states$values == island_state)
}

n_islands <- count_islands(decoded)
n_islands
