# paired end path:

deBruijnGraph_paired <- function(kmers) {
  # Initialize an empty list to store the adjacency list
  adjacencyList <- list()
  
  # Iterate through the kmers to extract paired prefixes and suffixes
  for (i in 1:length(kmers)) {
    pair <- strsplit(kmers[i], "\\|")[[1]]  # Split paired kmer by "|"
    k <- nchar(pair[1])  # Assuming both kmers in the pair have the same length
    
    prefix <- paste(substr(pair[1], 1, k - 1), substr(pair[2], 1, k - 1), sep=",")
    suffix <- paste(substr(pair[1], 2, k), substr(pair[2], 2, k), sep=",")
    
    # Update adjacency list treating prefix-suffix pairs as single entities
    if (!is.null(adjacencyList[[prefix]])) {
      adjacencyList[[prefix]] <- c(adjacencyList[[prefix]], suffix)
    } else {
      adjacencyList[[prefix]] <- c(suffix)
    }
  }
  
  # Return adjacency list
  adjacencyList
}


kmers<-c("ACC|ATA","ACT|ATT","ATA|TGA","ATT|TGA","CAC|GAT","CCG|TAC","CGA|ACT",
          "CTG|AGC","CTG|TTC","GAA|CTT","GAT|CTG","GAT|CTG","TAC|GAT","TCT|AAG","TGA|GCT","TGA|TCT","TTC|GAA")
# Generate overlap graph
result <- deBruijnGraph_paired(kmers)

eulerian_info <- find_eulerian_path_info(result)

# Execute Eulerian path computation
paths <- find_eulerian_paths(result, eulerian_info$Path_Length, eulerian_info$Start_Node)

# Output found paths
print(paste(paths[[1]], collapse = " "))

# Function to find the longest overlap
find_overlap <- function(seq1, seq2) {
  max_overlap <- 0
  overlap_seq <- ""
  
  for (i in 1:nchar(seq1)) {
    suffix <- substr(seq1, i, nchar(seq1))
    prefix <- substr(seq2, 1, nchar(suffix))
    
    if (suffix == prefix) {
      max_overlap <- nchar(suffix)
      overlap_seq <- suffix
      break
    }
  }
  
  return(overlap_seq)
}

# Function to reconstruct the final sequence
reconstruct_sequence <- function(paths) {
  first_seq <- c()
  second_seq <- c()
  
  for (entry in paths[[1]]) {
    parts <- strsplit(entry, ",")[[1]]
    first_seq <- c(first_seq, parts[1])
    second_seq <- c(second_seq, parts[2])
  }
  
  seq1 <- paste0(first_seq[1], paste(sapply(first_seq[-1], function(x) substr(x, nchar(x), nchar(x))), collapse = ""))
  seq2 <- paste0(second_seq[1], paste(sapply(second_seq[-1], function(x) substr(x, nchar(x), nchar(x))), collapse = ""))
  
  # Find the actual overlapping region
  overlap <- find_overlap(seq1, seq2)
  
  if (nchar(overlap) > 0) {
    seq2_trimmed <- substr(seq2, nchar(overlap) + 1, nchar(seq2))
    final_sequence <- paste0(seq1, seq2_trimmed)
  } else {
    final_sequence <- paste(seq1, seq2, sep = " ")
  }
  
  return(final_sequence)
}

reconstructed <- reconstruct_sequence(paths)
print(reconstructed)
