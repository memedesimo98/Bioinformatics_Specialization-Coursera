### peptide vector simplified

peptide <- "GQRFNILWYTFSYEINQVDFRKFHGHRLVGKSN"

Amino_map <- list(X = list(Mass = 4),Z = list(Mass = 5))

peptide_vectoring <- function(peptide,AA_map) {
  
  peptide_vector <- unlist(strsplit(peptide,""))
  
  peptide_vector_final <- c()
  
  for(i in 1:length(peptide_vector)) {
    amino <- peptide_vector[[i]]
    weight <- AA_map[[amino]]$Mass
    zeros <- rep(0,(weight - 1))
    zeros <- c(zeros,1)
    peptide_vector_final <- c(peptide_vector_final,zeros)
  }
  return(peptide_vector_final)
}

vector <- peptide_vectoring(peptide,AA_map)

vector <- unlist(strsplit(vector," "))

vector_to_peptide <- function(vector,AA_map) {
  
  if (is.character(vector)) {
    message("Input is a character vector")
    vector <- unlist(strsplit(vector," "))
  } else if (is.numeric(vector)) {
    # numeric covers both integer and double
    message("Input is numeric (integer or double)")
  }

  current_peptide <- c()
  peptide_vector <- character()
  
  while(length(vector)>0) {
    current_value <- vector[1]
    current_peptide <- c(current_peptide,current_value)
    vector <- vector[-1]
    if (current_value == 1) {
      actual_amino_length <- length(current_peptide)
      actual_amino <- names(AA_map)[
        sapply(AA_map, function(x) x$Mass) == actual_amino_length
      ][1]
      peptide_vector <- paste0(peptide_vector,actual_amino,collapse = "")
      current_peptide <- c()
    }
  }
  return(peptide_vector)
}

peptide <- vector_to_peptide(vector,AA_map)


# Amino acid map (toy alphabet)
AA_map <- list(X = list(Mass = 4), Z = list(Mass = 5))

# Function to score a peptide against a spectrum
score_peptide <- function(peptide, spectrum, AA_map) {
  peptide_vec <- peptide_vectoring(peptide, AA_map)
  if (length(peptide_vec) != length(spectrum)) {
    return(-Inf)  # discard peptides that don't match spectrum length
  }
  return(sum(peptide_vec * spectrum))
}


# Peptide Identification Problem solver
# Compute peptide length bounds based on spectrum and AA masses
length_bounds <- function(spectrum, AA_map) {
  L <- length(spectrum)
  masses <- sapply(AA_map, function(x) x$Mass)
  min_len <- ceiling(L / max(masses))
  max_len <- floor(L / min(masses))
  return(list(min_len = min_len, max_len = max_len))
}

# Peptide Identification with bounded loops
peptide_identification <- function(spectrum, proteome, AA_map) {
  bounds <- length_bounds(spectrum, AA_map)
  min_len <- bounds$min_len
  max_len <- bounds$max_len
  
  best_score <- -Inf
  best_peptide <- ""
  best_peptides <- character()
  
  for (i in 1:nchar(proteome)) {
    for (j in (i + min_len - 1):min(i + max_len - 1, nchar(proteome))) {
      candidate <- substr(proteome, i, j)
      sc <- score_peptide(candidate, spectrum, AA_map)
      
      if (sc > best_score) {
        best_score <- sc
        best_peptide <- candidate
        best_peptides <- candidate   # reset list
      } else if (sc == best_score && best_score != -Inf) {
        best_peptides <- c(best_peptides, candidate)
      }
    }
  }
  return(list(best_peptide = best_peptide, best_score = best_score, best_peptides = best_peptides))
}


# Function to read spectrum + proteome from file
read_spectrum_proteome <- function(filename) {
  # Read all lines
  lines <- readLines(filename)
  
  if (length(lines) < 2) {
    stop("File must contain at least two lines: spectrum and proteome")
  }
  
  # First line = spectrum (space-delimited numbers)
  spectrum <- as.numeric(unlist(strsplit(lines[1], " ")))
  
  # Second line = proteome string
  proteome <- lines[2]
  
  return(list(spectrum = spectrum, proteome = proteome))
}

# ---- Example usage ----
input_file <- choose.files()
data <- read_spectrum_proteome(input_file)
spectrum <- data$spectrum
proteome <- data$proteome
result <- peptide_identification(spectrum, proteome, AA_map)
print(result)


### fuller spectras

### test

# --- Utility functions ---

# Compute peptide length bounds based on spectrum and AA masses
length_bounds <- function(spectrum, AA_map) {
  L <- length(spectrum)
  masses <- sapply(AA_map, function(x) x$Mass)
  min_len <- ceiling(L / max(masses))
  max_len <- floor(L / min(masses))
  return(list(min_len = min_len, max_len = max_len))
}

# Convert peptide string into spectral vector
peptide_vectoring <- function(peptide, AA_map) {
  peptide_vector <- unlist(strsplit(peptide,""))
  peptide_vector_final <- c()
  
  for (i in seq_along(peptide_vector)) {
    amino <- peptide_vector[[i]]
    weight <- AA_map[[amino]]$Mass
    zeros <- rep(0, (weight - 1))
    zeros <- c(zeros, 1)
    peptide_vector_final <- c(peptide_vector_final, zeros)
  }
  return(peptide_vector_final)
}

# Score peptide against spectrum
score_peptide <- function(peptide, spectrum, AA_map) {
  peptide_vec <- peptide_vectoring(peptide, AA_map)
  if (length(peptide_vec) != length(spectrum)) {
    return(-Inf)  # discard peptides that don't match spectrum length
  }
  return(sum(peptide_vec * spectrum))
}

# Identify best peptide for a single spectrum
peptide_identification <- function(spectrum, proteome, AA_map) {
  bounds <- length_bounds(spectrum, AA_map)
  min_len <- bounds$min_len
  max_len <- bounds$max_len
  
  best_score <- -Inf
  best_peptide <- ""
  
  for (i in 1:nchar(proteome)) {
    for (j in (i + min_len - 1):min(i + max_len - 1, nchar(proteome))) {
      candidate <- substr(proteome, i, j)
      sc <- score_peptide(candidate, spectrum, AA_map)
      if (sc > best_score || sc == best_score) {
        best_score <- sc
        best_peptide <- candidate  # overwrite so "last wins" in ties
      }
    }
  }
  return(list(best_peptide = best_peptide, best_score = best_score))
}

library(parallel)

peptide_identification_over_spectra <- function(spectra, proteome, AA_map, threshold) {
  # Create a cluster using all available cores
  cl <- makeCluster(detectCores())
  
  # Export needed objects/functions to the cluster
  clusterExport(cl, c("proteome", "AA_map", "threshold",
                      "peptide_identification", "score_peptide",
                      "peptide_vectoring", "length_bounds"))
  
  # Run in parallel across spectra
  results <- parLapply(cl, spectra, function(current_spectrum) {
    current_result <- peptide_identification(current_spectrum, proteome, AA_map)
    if (current_result$best_score >= threshold) {
      return(current_result$best_peptide)
    } else {
      return(NULL)
    }
  })
  
  # Stop cluster
  stopCluster(cl)
  
  # Flatten results
  unlist(results)
}

input_file <- choose.files()
data <- read_spectrum_proteome(input_file)
spectra <- data$spectra
proteome <- data$proteome
threshold <- data$threshold
system.time({
  # Code you want to measure
  result <- peptide_identification_over_spectra(spectra, proteome, AA_map, threshold)
})
cat(result, sep = " ")
