#### viterbi most likely string (HMM)

most_likely_string <- function(n, alphabet, states, transition_mat, emission_mat) {
  k <- length(states)
  state_idx <- setNames(seq_len(k), states)
  
  # DP tables
  dp <- matrix(-Inf, nrow = n, ncol = k)
  bp_state <- matrix(NA_integer_, nrow = n, ncol = k)
  bp_symbol <- matrix(NA_character_, nrow = n, ncol = k)
  
  # Initialization: equal prior over states
  init <- rep(1/k, k)
  for (s in states) {
    si <- state_idx[[s]]
    # choose best emission symbol for this state
    best_sym <- alphabet[which.max(emission_mat[s, ])]
    dp[1, si] <- init[si] * max(emission_mat[s, ])
    bp_symbol[1, si] <- best_sym
  }
  
  # Recurrence
  for (i in 2:n) {
    for (s in states) {
      si <- state_idx[[s]]
      # best emission symbol for this state
      best_sym <- alphabet[which.max(emission_mat[s, ])]
      best_em <- max(emission_mat[s, ])
      
      candidates <- dp[i-1, ] * transition_mat[, s] * best_em
      best_prev <- which.max(candidates)
      
      dp[i, si] <- candidates[best_prev]
      bp_state[i, si] <- best_prev
      bp_symbol[i, si] <- best_sym
    }
  }
  
  # Termination
  last_best <- which.max(dp[n, ])
  
  # Traceback
  best_states <- character(n)
  best_symbols <- character(n)
  best_states[n] <- states[last_best]
  best_symbols[n] <- bp_symbol[n, last_best]
  
  for (i in (n-1):1) {
    prev <- bp_state[i+1, state_idx[best_states[i+1]]]
    best_states[i] <- states[prev]
    best_symbols[i] <- bp_symbol[i, prev]
  }
  
  list(
    emitted_string = paste(best_symbols, collapse=""),
    state_path = paste(best_states, collapse=""),
    prob = max(dp[n, ])
  )
}
