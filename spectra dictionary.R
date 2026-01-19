# spectra dictionary

size_spectral_dictionary <- function(spectrum, threshold, max_score, AA_map) {
  m <- length(spectrum)
  Size <- matrix(0, nrow = m+1, ncol = max_score+1)
  Size[1,1] <- 1  # empty peptide, score 0
  
  for (i in 1:m) {
    for (aa in names(AA_map)) {
      mass <- AA_map[[aa]]$Mass
      prev_i <- i - mass
      if (prev_i >= 0) {
        for (t in 0:max_score) {
          new_t <- t + spectrum[i]
          if (new_t >= 0 && new_t <= max_score) {
            Size[i+1, new_t+1] <- Size[i+1, new_t+1] + Size[prev_i+1, t+1]
          }
        }
      }
    }
  }
  
  return(sum(Size[m+1, (threshold+1):(max_score+1)]))
}



spec_str <- "13 -4 14 -2 1 5 15 6 -6 12 -9 14 10 -9 11 6 10 12 -2 15 5 8 -2 -4 -2 -7 14 -7 15 -1 9 -2 -8 0 -10 2 -9 3 -5 -9 -7 12 -7 0 10 4 -5 -6 15 11 10 6 1 7 -1 2 11 -2 -4 -10 -4 -4 -7 5 8 13 13 -6 -4 8 9 13 8 3 10 1 -6 13 3 0 1 -7 -9 -6 7 7 1 -9 -1 -7 11 9 4 11 13 -2 -2 10 15 -4 -6 11 -5 0 -3 -2 -4 -2 -9 -1 7 3 1 10 -4 14 1 -6 15 7 15 8 12 10 7 -10 8 3 8 -2 13 2 8 -8 1 9 2 -6 12 -1 14 15 -7 5 1 10 -6 -4 -7 1 0 10 15 14 5 -6 13 5 10 5 8 -3 2 10 -8 14 -4 -3 2 -2 12 9 5 11 1 -9 2 4 14 3 -1 -5 2 -8 -1 0 3 13 11 3 -6 0 12 -9 7 -8 5 -3 -8 2 -9 7 2 14 14 10 9 -5 -10 -8 0 15 -1 12 14 9 9 5 12 -8 5 -1 7 -8 -2 -8 1 9 -5 11 1 -8 -4 4 11 14 12 -8 -6 -1 -2 -5 12 -7 15 -6 -4 5 -4 9 14 1 -1 -4 8 13 3 13 -2 14 11 1 12 13 11 2 2 -8 -1 4 -7 -3 13 12 -4 -4 3 4 -9 4 13 -3 1 10 -4 14 -4 -9 11 -3 11 4 1 11 2 11 -4 15 -5 7 7 4 6 10 -3 3 5 8 15 6 0 13 15 -7 12 9 -1 -8 -2 11 15 10 -4 -8 2 8 -8 6 -1 -8 -5 14 -7 6 -10 15 -8 10 0 -5 2 -4 12 4 7 -9 -3 -8 8 -6 11 13 6 5 14 11 6 -8 15 8 6 -9 3 10 13 -1 12 -1 -8 -10 -3 9 11 6 1 9 -9 -5 15 4 4 7 12 -1 0 -4 -10 11 -10 -6 3 -10 8 9 -2 14 5 7 9 7 6 6 12 -10 -9 7 -9 -1 0 -8 7 -1 2 8 15 -2 8 -5 2 -2 3 3 15 14 1 -5 -4 -1 -9 12 13 8 -3 4 -8 4 -5 8 15 -6 4 -8 9 4 1 14 9 2 6 -1 -7 2 -3 -5 0 -8 1 -6 -1 -5 -6 12 -10 -6 15 -8 -4 -7 10 15 7 0 9 -5 13 2 13 15 0 3 7 7 7 -4 3 14 13 -5 13 11 1 9 -2 -10 -2 1 12 7 11 -2 5 -8 -7 7 -4 4 -4 6 1 13 11 -4 4 -8 1 -8 -2 4 -6 4 2 -3 13 8 -9"
spectrum <- as.integer(strsplit(spec_str, " ")[[1]])
threshold <- 35
max_score <- 200
result <- size_spectral_dictionary(spectrum, threshold, max_score, AA_map)
result


size_and_prob_spectral_dictionary <- function(spectrum, threshold, max_score, AA_map) {
  m <- length(spectrum)
  aa_masses <- vapply(AA_map, function(x) x$Mass, integer(1))
  p_aa <- 1.0 / length(aa_masses)  # uniform probability for each amino acid
  
  # DP tables: one for counts, one for probability mass
  Size <- matrix(0L, nrow = m+1, ncol = max_score+1)
  Prob <- matrix(0.0, nrow = m+1, ncol = max_score+1)
  
  # Base case: empty peptide at mass 0, score 0
  Size[1,1] <- 1L
  Prob[1,1] <- 1.0
  
  for (i in 1:m) {
    si <- spectrum[i]
    for (mass in aa_masses) {
      prev_i <- i - mass
      if (prev_i >= 0) {
        for (t in 0:max_score) {
          new_t <- t + si
          if (new_t >= 0 && new_t <= max_score) {
            # Update counts
            Size[i+1, new_t+1] <- Size[i+1, new_t+1] + Size[prev_i+1, t+1]
            # Update probability mass
            Prob[i+1, new_t+1] <- Prob[i+1, new_t+1] + Prob[prev_i+1, t+1] * p_aa
          }
        }
      }
    }
  }
  
  # Final results: counts and probability
  count <- sum(Size[m+1, (threshold+1):(max_score+1)])
  prob  <- sum(Prob[m+1, (threshold+1):(max_score+1)])
  
  return(list(count = count, probability = prob))
}

spec_str <- "7 -3 -2 -10 -4 2 8 8 10 2 -10 14 1 13 -9 3 2 2 15 -5 3 12 4 8 8 7 -5 -6 15 -5 14 -7 15 4 5 -6 -6 -2 12 -9 7 1 9 9 8 14 -1 8 -8 10 2 -4 -1 -7 2 -3 11 -5 -3 -2 -4 -7 -10 7 2 -9 -1 -6 0 -1 9 -9 14 11 9 -2 -4 4 -4 1 2 -2 7 12 2 13 -9 3 4 -3 1 3 -1 -2 14 1 4 0 7 11 7 10 -1 9 -9 -6 -2 -10 14 12 2 -6 4 0 -4 -8 -1 1 12 4 -5 -6 -9 -4 13 14 13 9 14 8 -9 -3 -1 4 6 3 -8 -2 -9 7 -8 -9 1 11 10 -10 11 -8 10 14 15 7 11 11 -3 -3 14 3 -8 -2 7 8 -5 2 8 1 2 9 8 -3 -1 -5 -8 15 1 15 6 13 10 4 -8 3 15 6 2 -10 -10 13 -10 -6 12 8 -8 13 -8 15 1 -9 14 15 11 -7 15 6 -7 6 -6 7 -3 -4 5 -5 -1 0 -1 8 8 -3 -6 -7 -1 -1 -3 -4 -4 6 -6 -3 10 8 7 5 -6 1 6 1 5 -4 -9 7 -10 7 0 7 -9 -4 -8 0 6 13 11 -9 -5 6 -5 11 -9 -4 0 14 11 7 8 -5 3 7 13 1 2 -5 2 4 -4 -6 8 1 -10 9 10 -9 8 0 -2 0 1 2 3 3 6 -6 3 3 5 14 0 -7 1 -2 12 -3 4 0 7 1 9 8 -3 4 -2 -3 11 4 -5 9 -8 5 -3 -7 -5 7 -8 0 -2 -3 -6 10 5 -10 3 13 -7 -1 12 4 -7 13 10 12 13 -2 12 9 14 7 6 -4 3 -10 -5 -10 1 3 4 3 13 14 11 15 5 8 13 -1 5 -2 -4 9 -1 6 -5 11 13 15 10 0 -3 -9 9 -7 14 1 -7 15 -2 8 1 -8 5 11 4 13 13 12 0 14 -4 1 -10 14 11 -5 3 -8 3 -4 -4 8 -6 -7 10 -4 -9 -9 -6 -3 15 11 -6 -9 10 -10 15 6 -9 3 -4 0 -8 -7 -8 3 -4 -9 8 -7 -4 11 15 -2 15 -9 2 0 -4 4 13 -2 -6 1 1 -5 -9 -6 2 -10 -10 1 10 14 0 4 15 -8 4 4 -7 -2 -1 -10 8 -5 -2 4 13 1 3 9 6 -2 5 11 9 6 10 15 8 6 8 14 6 0 13 2 13 -7 11 1 12 14 15 -6 12 12 -5 2 6 9 4 -6 4 3 3 5 -8 5 11 9 9 3 13 -10 7 -2"
spectrum <- as.integer(strsplit(spec_str, " ")[[1]])
threshold <- 34
max_score <- 200
result <- size_and_prob_spectral_dictionary(spectrum, threshold, max_score, AA_map)
result
