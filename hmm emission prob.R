### hmm emission prob

compute_emission_prob <- function(string, path, alphabet, states, emission_mat) {
  # Split into characters
  x <- strsplit(string, "")[[1]]
  pi <- strsplit(path, "")[[1]]
  
  if (length(x) != length(pi)) stop("String and path must have same length")
  
  prob <- 1
  for (i in seq_along(x)) {
    sym <- x[i]
    st  <- pi[i]
    prob <- prob * emission_mat[st, sym]
  }
  prob
}

# Example with your sample input
string <- "zzzyxyyzzx"
path <- "BAAAAAAAAA"
alphabet <- c("x","y","z")
states <- c("A","B")

emission_mat <- matrix(c(0.176,0.596,0.228,
                         0.225,0.572,0.203),
                       nrow=2, byrow=TRUE,
                       dimnames=list(states, alphabet))

prob <- compute_emission_prob(string, path, alphabet, states, emission_mat)
format(prob, digits=15)
# [1] "3.59748954746e-06"

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

# 3) Hidden path (line after second separator)
path <- trimws(lines[sep_idx[2] + 1])

# 4) States (line after third separator)
states <- unlist(strsplit(trimws(lines[sep_idx[3] + 1]), "\\s+"))

# 5) Emission matrix (block after fourth separator)
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

# Show results
cat("String:\n", string, "\n\n")
cat("Alphabet:\n", paste(alphabet, collapse = " "), "\n\n")
cat("Hidden path:\n", path, "\n\n")
cat("States:\n", paste(states, collapse = " "), "\n\n")
cat("Emission matrix:\n")
print(emission_mat)
