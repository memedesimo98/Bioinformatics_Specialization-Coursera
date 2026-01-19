### HMM viterbi forward probability

forward_prob <- function(x, alphabet, states, transition_mat, emission_mat) {
  x <- strsplit(x, "")[[1]]
  n <- length(x); k <- length(states)
  
  # Checks
  stopifnot(all(x %in% alphabet))
  
  # Equal initial probability
  init <- rep(1/k, k)
  
  # DP table: forward[i, s] = total prob of emitting prefix x1..xi ending in state s
  forward <- matrix(0, nrow = n, ncol = k)
  colnames(forward) <- states
  
  # Initialization
  sym1 <- x[1]
  for (s in states) {
    forward[1, s] <- init[states == s] * emission_mat[s, sym1]
  }
  
  # Recurrence
  for (i in 2:n) {
    sym <- x[i]
    for (s in states) {
      forward[i, s] <- sum(forward[i - 1, ] * transition_mat[, s]) * emission_mat[s, sym]
    }
  }
  
  # Termination: sum over ending states
  prob <- sum(forward[n, ])
  prob
}

x_str <- "TTHHH"
alphabet <- c("T","H")
states <- c("F","B")

transition_mat <- matrix(
  c(0.9, 0.1,
    0.1, 0.9),
  nrow = 2, byrow = TRUE,
  dimnames = list(states, states)
)

emission_mat <- matrix(
  c(0.5, 0.5,
    0.75, 0.25),
  nrow = 2, byrow = TRUE,
  dimnames = list(states, alphabet)
)

forward_prob(x_str, alphabet, states, transition_mat, emission_mat)
# Expected: 1.1005510319694847e-06
