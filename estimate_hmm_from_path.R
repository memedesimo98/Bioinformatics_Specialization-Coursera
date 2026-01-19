### estimate_hmm_from_path

estimate_hmm_from_path <- function(x,
                                   alphabet,
                                   states,
                                   path,
                                   pseudocount = 0,
                                   allowed_transitions = NULL,
                                   emitting_states = NULL) {
  # x: observed string (e.g., "AEFDFDC")
  # alphabet: character vector of symbols
  # states: character vector of all states (e.g., c("S","I0","M1","D1",...,"E"))
  # path: character vector of states from S-exclusive to E-exclusive (same format returned by viterbi_profile)
  # pseudocount: nonnegative scalar; 0 means MLE counts, >0 applies smoothing
  # allowed_transitions: optional matrix (states x states) of {0,1}; if NULL, all transitions are allowed
  # emitting_states: optional character vector of emitting states; if NULL, inferred as Mk and Ik (incl. I0)
  
  # --- validate inputs
  if (!is.character(x) || length(x) != 1) stop("x must be a single character string")
  seq_chars <- strsplit(x, "")[[1]]
  Lx <- length(seq_chars)
  if (!is.character(alphabet) || length(alphabet) == 0) stop("alphabet must be a non-empty character vector")
  if (!is.character(states) || length(states) == 0) stop("states must be a non-empty character vector")
  if (!is.character(path)) stop("path must be a character vector of states")
  if (!("S" %in% states)) stop("states must include 'S'")
  if (!("E" %in% states)) stop("states must include 'E'")
  if (any(!path %in% states)) stop("path contains states not in 'states'")
  if (pseudocount < 0) stop("pseudocount must be >= 0")
  
  # --- infer emitting vs silent states if not provided
  if (is.null(emitting_states)) {
    # Convention: emitting are M{int} and I{int} (including I0); silent are S, E, and D{int}
    is_M <- grepl("^M\\d+$", states)
    is_I <- grepl("^I\\d*$", states)
    is_D <- grepl("^D\\d+$", states)
    is_silent <- (states %in% c("S","E")) | is_D
    emitting_states <- states[!is_silent & (is_M | is_I)]
  } else {
    if (any(!emitting_states %in% states)) stop("emitting_states contains unknown states")
  }
  silent_states <- setdiff(states, emitting_states)
  
  # --- set up count matrices
  Tcounts <- matrix(0, nrow = length(states), ncol = length(states),
                    dimnames = list(states, states))
  Ecounts <- matrix(0, nrow = length(states), ncol = length(alphabet),
                    dimnames = list(states, alphabet))
  
  # --- build the full state sequence including S at start and E at end
  full_path <- c("S", path, "E")
  
  # --- transitions: count adjacent pairs along the path
  for (k in seq_len(length(full_path) - 1)) {
    from <- full_path[k]
    to   <- full_path[k + 1]
    Tcounts[from, to] <- Tcounts[from, to] + 1
  }
  
  # --- emissions: walk the path and consume x only at emitting states
  emit_positions <- which(full_path %in% emitting_states)
  if (length(emit_positions) != Lx) {
    stop(sprintf("Emission mismatch: path has %d emitting positions but x has length %d",
                 length(emit_positions), Lx))
  }
  for (i in seq_along(emit_positions)) {
    st <- full_path[emit_positions[i]]
    ch <- seq_chars[i]
    if (!ch %in% alphabet) stop(sprintf("Symbol '%s' not in alphabet", ch))
    Ecounts[st, ch] <- Ecounts[st, ch] + 1
  }
  
  # --- convert counts to probabilities with optional pseudocounts
  # Allowed transitions mask: if NULL, allow all transitions
  if (is.null(allowed_transitions)) {
    allowed_transitions <- matrix(1, nrow = length(states), ncol = length(states),
                                  dimnames = list(states, states))
  } else {
    # validate shape and names
    if (!all(dim(allowed_transitions) == c(length(states), length(states)))) {
      stop("allowed_transitions must have dim length(states) x length(states)")
    }
    if (!identical(rownames(allowed_transitions), states) ||
        !identical(colnames(allowed_transitions), states)) {
      stop("allowed_transitions dimnames must match 'states'")
    }
    # coerce to {0,1}
    allowed_transitions <- (allowed_transitions != 0) * 1
  }
  
  # transitions: add pseudocount only where allowed, then normalize each row
  Tmat <- Tcounts + pseudocount * allowed_transitions
  row_sums_T <- rowSums(Tmat)
  # If a row is completely forbidden (sum 0), leave as 0; else normalize
  nonzero_T <- row_sums_T > 0
  if (any(nonzero_T)) Tmat[nonzero_T, ] <- sweep(Tmat[nonzero_T, , drop = FALSE], 1, row_sums_T[nonzero_T], "/")
  if (any(!nonzero_T)) Tmat[!nonzero_T, ] <- 0
  
  # emissions: add pseudocount only on emitting rows; silent rows remain zero
  Emat <- Ecounts
  for (s in states) {
    if (s %in% emitting_states) {
      Emat[s, ] <- Emat[s, ] + pseudocount
      rs <- sum(Emat[s, ])
      if (rs > 0) {
        Emat[s, ] <- Emat[s, ] / rs
      } else {
        Emat[s, ] <- 0
      }
    } else {
      Emat[s, ] <- 0
    }
  }
  
  list(transition = Tmat, emission = Emat,
       counts = list(transition = Tcounts, emission = Ecounts),
       meta = list(emitting_states = emitting_states, silent_states = silent_states))
}

read_profile_input <- function(path) {
  lines <- trimws(readLines(path))
  # first non-empty line is threshold
  lines <- lines[nzchar(lines)]
  x <- lines[[1]]
  sep_idx <- which(lines == "--------")
  alphabet <- unlist(strsplit(lines[sep_idx[1] + 1]," "))
  path <- unlist(strsplit(lines[sep_idx[2] + 1], ""))
  states <- unlist(strsplit(lines[sep_idx[3] + 1], " "))
  states <- c("S",states,"E")
  list(x = x, alphabet = alphabet, path = path, states = states)
}

input_file <- choose.files()
input_lines <- read_profile_input(input_file)
hmm <- estimate_hmm_from_path(input_lines$x,input_lines$alphabet,input_lines$states,
                         input_lines$path)



### simpler test

read_profile_input <- function(path) {
  lines <- trimws(readLines(path))
  # first non-empty line is threshold
  lines <- lines[nzchar(lines)]
  x <- lines[[1]]
  sep_idx <- which(lines == "--------")
  alphabet <- unlist(strsplit(lines[sep_idx[1] + 1]," "))
  path <- unlist(strsplit(lines[sep_idx[2] + 1], ""))
  states <- unlist(strsplit(lines[sep_idx[3] + 1], " "))
  list(x = x, alphabet = alphabet, path = path, states = states)
}

estimate_simple_hmm <- function(x, alphabet, states, path) {
  seq_chars <- strsplit(x, "")[[1]]
  L <- length(seq_chars)
  if (length(path) != L) stop("Path length must equal sequence length")
  
  # Transition counts
  Tcounts <- matrix(0, nrow=length(states), ncol=length(states),
                    dimnames=list(states, states))
  for (i in seq_len(L-1)) {
    from <- path[i]
    to   <- path[i+1]
    Tcounts[from, to] <- Tcounts[from, to] + 1
  }
  
  # Emission counts
  Ecounts <- matrix(0, nrow=length(states), ncol=length(alphabet),
                    dimnames=list(states, alphabet))
  for (i in seq_len(L)) {
    st <- path[i]
    sym <- seq_chars[i]
    Ecounts[st, sym] <- Ecounts[st, sym] + 1
  }
  
  # Normalize
  Tmat <- Tcounts
  row_sums_T <- rowSums(Tmat)
  Tmat[row_sums_T > 0, ] <- sweep(Tmat[row_sums_T > 0, , drop=FALSE], 1, row_sums_T[row_sums_T > 0], "/")
  
  Emat <- Ecounts
  row_sums_E <- rowSums(Emat)
  Emat[row_sums_E > 0, ] <- sweep(Emat[row_sums_E > 0, , drop=FALSE], 1, row_sums_E[row_sums_E > 0], "/")
  
  list(transition=Tmat, emission=Emat)
}


input_file <- choose.files()
input_lines <- read_profile_input(input_file)
hmm <- estimate_simple_hmm(input_lines$x,input_lines$alphabet,input_lines$states,
                              input_lines$path)

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

