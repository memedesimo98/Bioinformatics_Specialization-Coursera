
viterbi_decode <- function(x, alphabet, states, transition_mat, emission_mat, start_prob = NULL) {
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
  
  # --- Start probability logic ---
  if (is.null(start_prob)) {
    # default: uniform
    start_prob <- rep(1/k, k)
    names(start_prob) <- states
  } else {
    stopifnot(all(names(start_prob) == states))
  }
  logInit <- log(start_prob)
  
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
    prob = exp(best_log_prob),
    dp = dp,
    bp = bp
  )
}

# Stati
states <- c("A+", "T+", "C+", "G+", "A-", "T-", "G-", "C-")

# Start probability (uniforme)
start_probability <- setNames(rep(0.125, length(states)), states)

# Transition probability matrix (p = 0.1)
p <- 0.1
transition_probability <- matrix(c(
  # A+
  0.180*(1-p), 0.122*(1-p), 0.268*(1-p), 0.430*(1-p), rep(0.25*p,4),
  # T+
  0.082*(1-p), 0.170*(1-p), 0.357*(1-p), 0.391*(1-p), rep(0.25*p,4),
  # C+
  0.191*(1-p), 0.211*(1-p), 0.299*(1-p), 0.299*(1-p), rep(0.25*p,4),
  # G+
  0.161*(1-p), 0.120*(1-p), 0.346*(1-p), 0.373*(1-p), rep(0.25*p,4),
  # A-
  rep(0.25*p,4), 0.300*(1-p), 0.210*(1-p), 0.200*(1-p), 0.290*(1-p),
  # T-
  rep(0.25*p,4), 0.176*(1-p), 0.291*(1-p), 0.242*(1-p), 0.291*(1-p),
  # C-
  rep(0.25*p,4), 0.319*(1-p), 0.291*(1-p), 0.302*(1-p), 0.081*(1-p),
  # G-
  rep(0.25*p,4), 0.251*(1-p), 0.199*(1-p), 0.251*(1-p), 0.299*(1-p)
), nrow=8, byrow=TRUE)

rownames(transition_probability) <- states
colnames(transition_probability) <- states

# Emission probability matrix
emission_probability <- matrix(0, nrow=8, ncol=4,
                               dimnames=list(states, c("A","T","G","C")))
emission_probability["A+", "A"] <- 1
emission_probability["T+", "T"] <- 1
emission_probability["G+", "G"] <- 1
emission_probability["C+", "C"] <- 1
emission_probability["A-", "A"] <- 1
emission_probability["T-", "T"] <- 1
emission_probability["G-", "G"] <- 1
emission_probability["C-", "C"] <- 1

# --- Stampa ---
x <- "GATCTGATAAGTCCCAGGACTTCAGAAGAGCTGTGAGACCTTGGCCAAGTCACTTCCTCCTTCAGGAACATTGCAGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAACATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGAAACTAAAGGGTCGCAGGCCGGTTTCTGCTAATTTCTTTAATTCCAAGACAGTCTCAAATATTTTCTTATTAACTTCCTGGAGGGAGGCTTATCATTCTCTCTTTTGGATGATTCTAAGTACCAGCTAAAATACAGCTATCATTCATTTTCCTTGATTTGGGAGCCTAATTTCTTTAATTTAGTATGCAAGAAAACCAATTTGGAAATATCAACTGTTTTGGAAACCTTAGACCTAGGTCATCCTTAGTAAGATCTTCCCATTTATATAAATACTTGCAAGTAGTAGTGCCATAATTACCAAACATAAAGCCAACTGAGATGCCCAAAGGGGGCCACTCTCCTTGCTTTTCCTCCTTTTTAGAGGATTTATTTCCCATTTTTCTTAAAAAGGAAGAACAAACTGTGCCCTAGGGTTTACTGTGTCAGAACAGAGTGTGCCGATTGTGGTCAGGACTCCATAGCATTTCACCATTGAGTTATTTCCGCCCCCTTACGTGTCTCTCTTCAGCGGTCTATTATCTCCAAGAGGGCATAAAACACTGAGTAAACAGCTCTTTTATATGTGTTTCCTGGATGAGCCTTCTTTTAATTAATTTTGTTAAGGGATTTCCTCTAGGGCCACTGCACGTCATGGGGAGTCACCCCCAGACACTCCCAATTGGCCCCTTGTCACCCAGGGGCACATTTCAGCTATTTGTAAAACCTGAAATCACTAGAAAGGAATGTCTAGTGACTTGTGGGGGCCAAGGCCCTTGTTATGGGGATGAAGGCTCTTAGGTGGTAGCCCTCCAAGAGAATAGATGGTG"
alphabet <- c("A","T","G","C")
result <- viterbi_decode(x,alphabet,states,transition_probability,
                         emission_probability,start_probability)
x <- "ATGGCCCGAACCAAGCAGACTGCGCGCAAGTCAACGGGTGGCAAGGCGCCGCGCAAGCAGCTGGCCACCAAGGTGGCTCGCAAGAGCGCACCTGCCACTGGCGGCGTGAAGAAGCCGCACCGCTACCGGCCCGGCACGGTGGCGCTTCGCGAGATCCGCCGCTACCAGAAGTCCACTGAGCTGCTAATCCGCAAGTTGCCCTTCCAGCGGCTGATGCGCGAGATCGCTCAGGACTTTAAGACCGACCTGCGCTTCCAGAGCTCGGCCGTGATGGCGCTGCAGGAGGCGTGCGAGTCTTACCTGGTGGGGCTGTTTGAGGACACCAACCTGTGTGTCATCCATGCCAAACGGGTCACCATCATGCCTAAGGACATCCAGCTGGCACGCCGTATCCGCGGGGAGCGGGCCTAGGAGGGCTATCTCGCCACCTGAGAGGTTGCGCAACGTTCACCCCAAAGGCTCTTTTAAGAGCCACCCACCT"
x <- unlist(strsplit(x,""))
cat(paste(x,sep = ","))
