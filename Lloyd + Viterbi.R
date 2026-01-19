### Lloyd + Viterbi

read_profile_input <- function(path) {
  lines <- trimws(readLines(path))
  # first non-empty line is threshold
  lines <- lines[nzchar(lines)]
  iteration <- as.numeric(lines[1])
  # find separators '--------'
  sep_idx <- which(lines == "--------")
  x <- lines[sep_idx[1] + 1]
  alphabet <- strsplit(lines[sep_idx[2] + 1], "\\s+")[[1]]
  states <- strsplit(lines[sep_idx[3] + 1], "\\s+")[[1]]
  
  # 4) Emission matrix (block after third separator)
  mat_start <- sep_idx[4] + 1
  mat_end <- sep_idx[5] - 1
  mat_lines <- lines[mat_start:mat_end]
  
  transition_mat <- as.matrix(read.table(text = mat_lines, header = TRUE, row.names = 1))
  
  # 5) Emission matrix (block after 4th separator)
  mat_start <- sep_idx[5] + 1
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
  list(iteration = iteration, alphabet = alphabet, states = states, transition = transition_mat,
       x = x,emission = emission_mat)
}

simple_viterbi_from_hmm <- function(hmm, x) {
  states <- rownames(hmm$transition)
  trans <- hmm$transition
  emit  <- hmm$emission
  
  seq_chars <- strsplit(x, "")[[1]]
  L <- length(seq_chars)
  nstates <- length(states)
  
  # log probabilities (avoid underflow)
  log_trans <- log(trans)
  log_trans[!is.finite(log_trans)] <- -Inf
  log_emit <- log(emit)
  log_emit[!is.finite(log_emit)] <- -Inf
  
  # DP matrices
  V <- matrix(-Inf, nrow=nstates, ncol=L, dimnames=list(states, 1:L))
  back <- matrix(NA_character_, nrow=nstates, ncol=L, dimnames=list(states, 1:L))
  
  # initialization: first symbol
  ch <- seq_chars[1]
  for (s in states) {
    V[s,1] <- log_emit[s,ch]   # assume uniform start prob
  }
  
  # recursion
  for (t in 2:L) {
    ch <- seq_chars[t]
    for (s in states) {
      le <- log_emit[s,ch]
      if (is.infinite(le)) next
      best_val <- -Inf
      best_prev <- NA_character_
      for (p in states) {
        cand <- V[p,t-1] + log_trans[p,s] + le
        if (cand > best_val) {
          best_val <- cand
          best_prev <- p
        }
      }
      V[s,t] <- best_val
      back[s,t] <- best_prev
    }
  }
  
  # termination: best state at last position
  last_probs <- V[,L]
  cur_state <- states[which.max(last_probs)]
  
  # backtrack
  path <- character(L)
  path[L] <- cur_state
  for (t in (L-1):1) {
    cur_state <- back[cur_state,t+1]
    path[t] <- cur_state
  }
  
  return(path)
}


Viterbi_Lloyd <- function(transition,emission,x,alphabet,states,iteration = 100) {
  
  i <- 0
  
  while (i < iteration) {
    hmm <- list(transition = transition, emission = emission)
    path_identification <- simple_viterbi_from_hmm(hmm,x)
    Matrixes_creation <- estimate_simple_hmm(x,alphabet,states,path_identification)
    transition <- Matrixes_creation[[1]]
    emission <- Matrixes_creation[[2]]
    i <- i + 1
  }
    list(transition = transition, emission = emission)
}

input_file <- choose.files()
input_lines <- read_profile_input(input_file)
hmm <- Viterbi_Lloyd(input_lines$transition, input_lines$emission,
                     input_lines$x,input_lines$alphabet,input_lines$states,
                           input_lines$iteration)

con <- file("result.txt", "w", encoding = "UTF-8")
cat("\t", paste(colnames(hmm$transition), collapse = "\t"), "\n", file = con, sep = "")
for(r in rownames(hmm$transition)) {
  cat(r, "\t", paste(sprintf("%.3f", hmm$transition[r, ]), collapse = "\t"), "\n", file = con, sep = "")
}
cat("--------\n", file = con, sep = "")
cat("\t", paste(colnames(hmm$emission), collapse = "\t"), "\n", file = con, sep = "")
for(r in rownames(hmm$emission)) {
  cat(r, "\t", paste(sprintf("%.3f", hmm$emission[r, ]), collapse = "\t"), "\n", file = con, sep = "")
}
close(con)
