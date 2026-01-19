### hmm hidden path

compute_hidden_path_prob <- function(path, states, transition_mat) {
  # path: string like "ABABBBAAAA"
  # states: character vector of state names, e.g., c("A","B")
  # transition_mat: matrix with row/col names equal to states
  
  # Initial probability: uniform over states
  p <- 1 / length(states)
  
  # Iterate over consecutive pairs in path
  chars <- strsplit(path, "")[[1]]
  for (i in seq_len(length(chars) - 1)) {
    from <- chars[i]
    to   <- chars[i + 1]
    p <- p * transition_mat[from, to]
  }
  p
}

# Example using the sample input
states <- c("A", "B")
transition_mat <- matrix(c(0.377, 0.623,
                           0.26,  0.74),
                         nrow = 2, byrow = TRUE,
                         dimnames = list(states, states))
path <- "ABABBBAAAA"

prob <- compute_hidden_path_prob(path, states, transition_mat)
format(prob, digits = 15)
# [1] "0.000384928691755"


# Lazy-friendly HMM input loader
input_file <- choose.files()

# Read file lines
lines <- readLines(input_file)

# First line = hidden path
path <- lines[1]

# Second block = states (split by space)
states <- unlist(strsplit(lines[3], "\\s+"))

# Remaining lines = transition matrix
mat_lines <- lines[6:length(lines)]
transition_mat <- read.table(text = mat_lines, header = TRUE, row.names = 1)

# Check results
cat("Hidden path:\n", path, "\n\n")
cat("States:\n", states, "\n\n")
cat("Transition matrix:\n")
print(transition_mat)


### test

# Pick file
input_file <- choose.files()

# Read and clean lines
lines <- readLines(input_file, warn = FALSE)
lines <- lines[nchar(trimws(lines)) > 0]  # drop empty lines

# Locate separators like "--------"
is_sep <- grepl("^-{2,}$", trimws(lines))
sep_idx <- which(is_sep)
if (length(sep_idx) < 2) stop("Expected two '--------' separators.")

# 1) Hidden path (line before first separator)
path <- trimws(lines[1])

# 2) States (line after first separator)
states_line <- trimws(lines[sep_idx[1] + 1])
states <- unlist(strsplit(states_line, "\\s+"))

# After locating matrix block
mat_lines <- lines[mat_start:length(lines)]

# Drop the first line (column headers "A B")
mat_lines <- mat_lines[-1]

# Now parse
mat_df <- read.table(text = mat_lines,
                     header = FALSE, stringsAsFactors = FALSE,
                     sep = "", fill = TRUE, strip.white = TRUE)

rownames_vec <- mat_df[[1]]
numeric_mat <- as.matrix(mat_df[, -1, drop = FALSE])
storage.mode(numeric_mat) <- "numeric"
rownames(numeric_mat) <- rownames_vec
colnames(numeric_mat) <- states

transition_mat <- numeric_mat
print(transition_mat)

transition_mat <- numeric_mat

# Probability of hidden path with uniform initial state
compute_hidden_path_prob <- function(path, states, transition_mat) {
  chars <- strsplit(path, "")[[1]]
  if (!all(chars %in% states)) {
    stop("Path contains states not listed in 'states'.")
  }
  p <- 1 / length(states)
  for (i in seq_len(length(chars) - 1)) {
    from <- chars[i]
    to   <- chars[i + 1]
    p <- p * transition_mat[from, to]
  }
  p
}

# Run and print
prob <- compute_hidden_path_prob(path, states, transition_mat)

cat("Hidden path:\n", path, "\n\n")
cat("States:\n", paste(states, collapse = " "), "\n\n")
cat("Transition matrix:\n")
print(transition_mat)
cat("\nProbability of path:\n", format(prob, digits = 15), "\n")
