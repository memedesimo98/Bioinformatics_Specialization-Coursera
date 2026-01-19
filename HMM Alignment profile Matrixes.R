### HMM Alignment profile Matrixes

# helper to print matrices in the challenge format
print_matrix <- function(mat) {
  cat("\t", paste(colnames(mat), collapse = "\t"), "\n")
  for(r in rownames(mat)) {
    cat(r, "\t", paste(sprintf("%.3f", mat[r, ]), collapse = "\t"), "\n")
  }
}



read_profile_input <- function(path) {
  lines <- trimws(readLines(path))
  # first non-empty line is threshold
  lines <- lines[nzchar(lines)]
  theta <- as.numeric(unlist(strsplit(lines[1]," ")))
  pseudocounts <- theta[[2]]
  theta <- theta[[1]]
  # find separators '--------'
  sep_idx <- which(lines == "--------")
  if(length(sep_idx) < 2) stop("Unexpected input format: need two '--------' lines")
  alphabet <- strsplit(lines[sep_idx[1] + 1], "\\s+")[[1]]
  seq_profile <- lines[(sep_idx[2] + 1):length(lines)]
  seq_profile <- seq_profile[nzchar(seq_profile)]
  list(threshold = theta, alphabet = alphabet, seq_profile = seq_profile, pseudocounts = pseudocounts)
}

build_profile_hmm <- function(threshold, alphabet, seq_profile, pseudocount = 0.01) {
  seq_profile <- seq_profile[nzchar(seq_profile)]
  n <- length(seq_profile)
  if(n == 0) stop("Empty sequence profile")
  L <- nchar(seq_profile[1])
  aln <- do.call(rbind, lapply(seq_profile, function(s) strsplit(s, "")[[1]]))
  gap_frac <- apply(aln, 2, function(col) sum(col == "-") / n)
  is_match_col <- gap_frac < threshold
  M <- sum(is_match_col)
  
  # build state list
  states <- c("S", "I0")
  for(i in seq_len(M)) states <- c(states, paste0("M", i), paste0("D", i), paste0("I", i))
  states <- c(states, "E")
  
  # emitting states (M and I, including I0)
  emit_states <- states[grepl("^M\\d+$|^I\\d*$", states)]
  
  # initialize counts and allowed-transitions mask
  Tcounts <- matrix(0, nrow = length(states), ncol = length(states),
                    dimnames = list(states, states))
  Ecounts <- matrix(0, nrow = length(states), ncol = length(alphabet),
                    dimnames = list(states, alphabet))
  allowed <- matrix(0, nrow = length(states), ncol = length(states),
                    dimnames = list(states, states))
  
  # helper to mark allowed edges
  allow <- function(from, to) {
    if(from %in% rownames(allowed) && to %in% colnames(allowed)) allowed[from, to] <<- 1
  }
  
  # allowed edges according to profile HMM diagram
  allow("S", "I0")
  if(M >= 1) { allow("S", "M1"); allow("S", "D1") } else { allow("S", "E") }
  
  allow("I0", "I0")
  if(M >= 1) { allow("I0", "M1"); allow("I0", "D1") }
  
  for(k in seq_len(M)) {
    Mk <- paste0("M", k); Dk <- paste0("D", k); Ik <- paste0("I", k)
    # self-insertions
    allow(Ik, Ik)
    # transitions to insertion
    allow(Mk, Ik); allow(Dk, Ik)
    # transitions to next match/delete or End
    if(k < M) {
      nextM <- paste0("M", k+1); nextD <- paste0("D", k+1)
      allow(Mk, nextM); allow(Mk, nextD)
      allow(Dk, nextM); allow(Dk, nextD)
      allow(Ik, nextM); allow(Ik, nextD)
    } else {
      allow(Mk, "E"); allow(Dk, "E"); allow(Ik, "E")
    }
  }
  
  # count transitions and emissions from alignment (do not gate on allowed)
  for(r in seq_len(n)) {
    cur_state <- "S"
    match_idx <- 0
    for(c in seq_len(L)) {
      ch <- aln[r, c]
      if(is_match_col[c]) {
        match_idx <- match_idx + 1
        if(ch == "-") {
          next_state <- paste0("D", match_idx)
        } else {
          next_state <- paste0("M", match_idx)
          Ecounts[next_state, ch] <- Ecounts[next_state, ch] + 1
        }
      } else {
        if(ch == "-") {
          next_state <- NULL
        } else {
          next_state <- paste0("I", match_idx)
          Ecounts[next_state, ch] <- Ecounts[next_state, ch] + 1
        }
      }
      if(!is.null(next_state)) {
        Tcounts[cur_state, next_state] <- Tcounts[cur_state, next_state] + 1
        cur_state <- next_state
      }
    }
    # final transition to End
    Tcounts[cur_state, "E"] <- Tcounts[cur_state, "E"] + 1
  }
  
  # convert counts -> probability matrices (no pseudocount yet)
  Tmat <- Tcounts
  row_sums_T <- rowSums(Tmat)
  nonzero_T <- row_sums_T > 0
  if(any(nonzero_T)) Tmat[nonzero_T, ] <- Tmat[nonzero_T, ] / row_sums_T[nonzero_T]
  Tmat[!nonzero_T, ] <- 0
  
  Emat <- Ecounts
  row_sums_E <- rowSums(Emat)
  nonzero_E <- row_sums_E > 0
  if(any(nonzero_E)) Emat[nonzero_E, ] <- Emat[nonzero_E, ] / row_sums_E[nonzero_E]
  Emat[!nonzero_E, ] <- 0
  
  # add pseudocount sigma to allowed transition entries only (forbidden remain zero)
  Tmat <- Tmat + pseudocount * allowed
  # renormalize transition rows
  row_sums_T <- rowSums(Tmat)
  nonzero_T <- row_sums_T > 0
  if(any(nonzero_T)) Tmat[nonzero_T, ] <- Tmat[nonzero_T, ] / row_sums_T[nonzero_T]
  Tmat[!nonzero_T, ] <- 0
  
  # add pseudocount sigma to emission probabilities only for emitting states, then renormalize each row
  for(s in rownames(Emat)) {
    if(s %in% emit_states) {
      Emat[s, ] <- Emat[s, ] + pseudocount
      rs <- sum(Emat[s, ])
      if(rs > 0) Emat[s, ] <- Emat[s, ] / rs else Emat[s, ] <- 0
    } else {
      Emat[s, ] <- 0
    }
  }
  
  list(transition = Tmat, emission = Emat)
}


# Example: read from stdin
input_file <- choose.files()
input_lines <- read_profile_input(input_file)
hmm <- build_profile_hmm(input_lines$threshold,input_lines$alphabet,input_lines$seq_profile,
                         input_lines$pseudocounts)


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

