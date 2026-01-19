### spectra bulding improvements

peptide <- "MA"

# Choose the file and store the path
file_path <- file.choose("text")

# Open and read the content of the text file
file_content <- readLines(file_path)

practical_spectrum <- as.numeric(unlist(strsplit(file_content," ")))

peptide_spectrum <- Peptide_cyclic_Spectrum(peptide,AA_map)

practical_spectrum <- peptide_spectrum
score <- check_similarity_in_spectra(peptide_spectrum,practical_spectrum)

check_similarity_in_spectra <- function (peptide_spectrum,practical_spectrum) {
  score <- 0
  for (mass in peptide_spectrum) {
    if (mass %in% practical_spectrum) {
      set_of_masses <- practical_spectrum[practical_spectrum == mass]
      set_of_masses <- set_of_masses[-1]
      practical_spectrum <- practical_spectrum[practical_spectrum != mass]
      if (length(set_of_masses)>0) {
        practical_spectrum <- c(practical_spectrum,set_of_masses)
      }
      score <- score + 1
    }
  }
  return (score)
}

remove_cyclic_variants <- function(x) {
  unique_elements <- character(0)  # Initialize an empty vector
  
  for (element in x) {
    n <- nchar(element)  # Get length
    cycle_seq <- substr(paste0(element, element), 1, n + n - 1)  # Create cyclic version minus last char
    
    cyclic_variants <- unique(sapply(1:n, function(i) substr(cycle_seq, i, i + n - 1)))  # Extract cyclic versions

    # Check if any existing unique element matches a cyclic variant
    if (!any(unique_elements %in% cyclic_variants)) {
      unique_elements <- c(unique_elements, element)
    }
  }
  
  return(unique_elements)
}

subpeptides_scoring_finder <- function(peptide,N,AA_map) {
  if (is.character(peptide)) {
    peptide <- Peptide_cyclic_Spectrum(peptide,AA_map)
  }
  practical_spectrum <- peptide
  final_subpeptides <- character()
  candidates <- ""
  best_leaderboard_score <- - Inf
  
  while(length(candidates) > 0) {
    candidates <- sapply(candidates, function(x) {
      new_candidates <- character()
      new_candidates <- sapply(names(AA_map),function(aa) {
        new_candidates<-c(new_candidates,paste0(x,aa))
        })
      return(new_candidates)
    })
    
    candidates <- remove_cyclic_variants(candidates)
    
    print("candidates")
    print(candidates[1])
    # Generate candidate spectra as a named list
    candidates_spectra <- setNames(
      lapply(candidates, function(seq) Peptide_cyclic_Spectrum(seq, AA_map)),
      candidates
    )
    
    pre_final_candidates <- lapply(candidates_spectra, function(spectrum) 
      check_similarity_in_spectra(spectrum, practical_spectrum))
    
    leaderboard_pre_final_candidates <- rev(sort(sapply(pre_final_candidates,unique)))
    leaderboard_trim <- unique(leaderboard_pre_final_candidates[1:N])
    print("leaderboard_trim")
    print(leaderboard_trim)
    if (any(leaderboard_trim > best_leaderboard_score)) {
      best_leaderboard_score <- leaderboard_trim[1]
      print("best_leaderboard_score")
      print(best_leaderboard_score)
    final_candidates <- names(pre_final_candidates[pre_final_candidates %in% leaderboard_trim])
    final_subpeptides <- final_candidates
    
    candidates <- final_candidates
    } else {
    candidates <- character()
    }
  }
  return(final_subpeptides)
}

# Choose the file and store the path
file_path <- file.choose("text")

# Open and read the content of the text file
file_content <- readLines(file_path)
file_content <- "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
practical_spectrum <- as.numeric(unlist(strsplit(file_content," ")))

peptide <- practical_spectrum
n <- 151



og_time <- system.time(result <- subpeptides_scoring_finder(peptide,n,AA_map))
parallel_time <- system.time(result_parallel <- subpeptides_scoring_finder_parallel(peptide,n,AA_map))
result <- subpeptides_scoring_finder(peptide,n,AA_map)
final_masses <- subpeptide_to_masses(result,AA_map)
paste(unique(final_masses)[1],collapse = " ")

install.packages("foreach")
install.packages("doParallel")
library(foreach)
library(parallel)
library(doParallel)



subpeptides_scoring_finder_parallel <- function(peptide,N,AA_map) {
  if (is.character(peptide)) {
    peptide <- Peptide_cyclic_Spectrum(peptide,AA_map)
  }
  practical_spectrum <- peptide
  final_subpeptides <- character()
  candidates <- ""
  best_leaderboard_score <- - Inf
  # Define the number of cores
  num_cores <- 8
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary functions and variables just once
  clusterExport(cl, c("Peptide_cyclic_Spectrum", "check_similarity_in_spectra", "AA_map", "practical_spectrum"))
  
  # Start the while loop
  while (length(candidates) > 0) {
    
    candidates <- foreach(x = candidates, .combine = c) %dopar% {
      new_candidates <- unlist(lapply(names(AA_map), function(aa) paste0(x, aa)))
      return(new_candidates)
    }
    candidates <- remove_cyclic_variants(candidates)
    # Parallel processing for spectra generation
    candidates_spectra_list <- foreach(seq = candidates, .combine = "c") %dopar% {
      list(Peptide_cyclic_Spectrum(seq, AA_map))  # Ensure each output remains a list
    }
    candidates_spectra <- candidates_spectra_list
    # Parallel processing for similarity checking
    pre_final_candidates_list <- foreach(spectrum = candidates_spectra, .combine = c) %dopar% {
      check_similarity_in_spectra(spectrum, practical_spectrum)
    }
    names(pre_final_candidates_list) <- candidates
    # Process leaderboard logic
    leaderboard_pre_final_candidates <- rev(sort(sapply(pre_final_candidates_list, unique)))
    
    if (length(leaderboard_pre_final_candidates) < n) {
      leaderboard_trim <- unique(leaderboard_pre_final_candidates)
    } else {
    leaderboard_trim <- unique(leaderboard_pre_final_candidates[1:N])
    }
    if (any(leaderboard_trim > best_leaderboard_score)) {
      best_leaderboard_score <- leaderboard_trim[1]
      final_candidates <- names(pre_final_candidates_list[pre_final_candidates_list %in% leaderboard_trim])
      final_subpeptides <- final_candidates
      candidates <- final_candidates
    } else {
      candidates <- character()
    }
  }
  
  # Stop the cluster after finishing all iterations
  stopCluster(cl)
  
  return(final_subpeptides)
}

subpeptides_scoring_finder_parallel <- function(peptide,N,AA_map) {
  if (is.character(peptide)) {
    peptide <- Peptide_cyclic_Spectrum(peptide,AA_map)
  }
  practical_spectrum <- peptide
  final_subpeptides <- character()
  candidates <- ""
  best_leaderboard_score <- - Inf
  # Define the number of cores
  num_cores <- 8
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary functions and variables just once
  clusterExport(cl, c("Peptide_cyclic_Spectrum", "check_similarity_in_spectra", "AA_map", "practical_spectrum"))
  
  # Start the while loop
  while (length(candidates) > 0) {
    
    candidates <- foreach(x = candidates, .combine = c) %dopar% {
      new_candidates <- unlist(lapply(names(AA_map), function(aa) paste0(x, aa)))
      return(new_candidates)
    }
    candidates <- remove_cyclic_variants(candidates)
    # Parallel processing for spectra generation
    candidates_spectra_list <- foreach(seq = candidates, .combine = "c") %dopar% {
      list(Peptide_cyclic_Spectrum(seq, AA_map))  # Ensure each output remains a list
    }
    candidates_spectra <- candidates_spectra_list
    # Parallel processing for similarity checking
    pre_final_candidates_list <- foreach(spectrum = candidates_spectra, .combine = c) %dopar% {
      check_similarity_in_spectra(spectrum, practical_spectrum)
    }
    names(pre_final_candidates_list) <- candidates
    # Process leaderboard logic
    leaderboard_pre_final_candidates <- rev(sort(sapply(pre_final_candidates_list, unique)))
    
    if (length(leaderboard_pre_final_candidates) < n) {
      leaderboard_trim <- unique(leaderboard_pre_final_candidates)
    } else {
      leaderboard_trim <- unique(leaderboard_pre_final_candidates[1:N])
    }
    if (any(leaderboard_trim > best_leaderboard_score)) {
      best_leaderboard_score <- leaderboard_trim[1]
      final_candidates <- names(pre_final_candidates_list[pre_final_candidates_list %in% leaderboard_trim])
      final_subpeptides <- final_candidates
      candidates <- final_candidates
    } else {
      candidates <- character()
    }
  }
  
  # Stop the cluster after finishing all iterations
  stopCluster(cl)
  
  return(final_subpeptides)
}
# Choose the file and store the path
file_path <- file.choose("text")

# Open and read the content of the text file
file_content <- readLines(file_path)
practical_spectrum <- as.numeric(unlist(strsplit(file_content," ")))
peptide <- practical_spectrum
n <- 1000
parallel_time <- system.time(result_parallel <- subpeptides_scoring_finder_parallel(peptide,n,AA_map))
final_masses <- subpeptide_to_masses(result_parallel,AA_map)
paste(unique(final_masses)[1],collapse = " ")
parallel_time
