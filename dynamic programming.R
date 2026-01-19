# from m to n of peptides

ReduceAAmap <- function(AA_map) {
  filtered_map <- list()  # Store filtered entries
  seen_masses <- numeric(0)  # Track unique masses
  
  for (name in names(AA_map)) {
    mass <- AA_map[[name]]$Mass
    
    # Ensure valid mass before proceeding
    if (length(mass) == 0) next
    
    # Keep only the first occurrence of each unique mass
    if (!(mass %in% seen_masses)) {
      seen_masses <- c(seen_masses, mass)  # Add mass to tracked list
      filtered_map[[name]] <- AA_map[[name]]  # Keep entry
    }
  }
  
  return(filtered_map)
}

reducted_AA_map <- ReduceAAmap(AA_map)

trova_numero_combinazioni <- function(aminoacidi, M) {
  dp <- integer(M + 1)
  dp[1] <- 1  # dp[0] = 1
  
  for (massa in 0:M) {
    if (dp[massa + 1] > 0) {
      for (peso in aminoacidi) {
        nuova_massa <- massa + peso
        if (nuova_massa <= M) {
          dp[nuova_massa + 1] <- dp[nuova_massa + 1] + dp[massa + 1]
        }
      }
    }
  }
  
  return(dp[M + 1])
}

# Esempio di uso:

aminoacidi <- c(A = 2, B = 3, C = 5)
M <- 1024
named_mass_vector <- sapply(reducted_AA_map, function(x) x$Mass)
risultato <- trova_numero_combinazioni(named_mass_vector, M)
print(format(risultato,scientific = FALSE))




# calculate n* of subpeptides

n <- 34024

subpeptide_calculator <- function(n){
  subpeptides <- 1
  for (sub_n in 0:(n-1)) {
    subpeptides <- subpeptides + n - sub_n
  }
  return(subpeptides)
}

subpeptites <- subpeptide_calculator(n)
