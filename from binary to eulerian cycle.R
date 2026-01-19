deBruijnGraph_test <- function(text) {
  # Initialize an empty list to store the adjacency list
  adjacencyList <- list()
  
  # Iterate through the text to create k-mers and (k-1)-mers
  for (i in 1:length(text) ) {
    k<-nchar(kmers[i])
    kmer <- kmers[i]
    prefix <- substr(kmer, 1, k - 1)
    suffix <- substr(kmer, 2, k)
    
    # Update adjacency list
    if (!is.null(adjacencyList[[prefix]])) {
      adjacencyList[[prefix]] <- c(adjacencyList[[prefix]], suffix)
    } else {
      adjacencyList[[prefix]] <- c(suffix)
    }
  }
  
  # Return adjacency list
  adjacencyList
}

# Function to generate all binary strings of length k
generate_binary_combinations <- function(k) {
  # Create all possible combinations using expand.grid()
  binary_matrix <- expand.grid(rep(list(c("0", "1")), k))
  
  # Convert matrix rows into character strings
  binary_strings <- apply(binary_matrix, 1, paste0, collapse = "")
  
  return(binary_strings)
}

find_eulerian_cycle <- function(data) {
  # Parse the input data to build the graph
  graph <- list()
  for (entry in data) {
    parts <- strsplit(entry, ": ")[[1]]
    node <- parts[1]  # Keep node as a character
    if (length(parts) > 1) {
      neighbors <- unlist(strsplit(parts[2], " "))
    } else {
      neighbors <- character(0)
    }
    graph[[node]] <- neighbors
  }
  
  # Initialize variables
  cycle <- character(0)
  stack <- c(names(graph)[1]) # Start at any node
  
  # Traverse the graph using Hierholzer's algorithm
  while (length(stack) > 0) {
    current <- stack[length(stack)]
    if (length(graph[[current]]) > 0) {
      # Follow an edge
      next_node <- graph[[current]][1]
      graph[[current]] <- graph[[current]][-1]  # Remove used edge
      stack <- c(stack, next_node)
    } else {
      # No more neighbors; add to cycle
      cycle <- c(cycle, current)
      stack <- stack[-length(stack)]
    }
  }
  
  # Reverse the cycle to get the correct order
  return(rev(cycle))  # Ensure we return character output
}

# Function to reconstruct the original sequence from edge order
reconstruct_string <- function(edge_order) {
  # Start with the first node
  reconstructed_string <- edge_order[1]
  k<-nchar(reconstructed_string)
  # Iterate through the remaining nodes
  for (i in 2:length(edge_order)) {
    next_node <- edge_order[i]
    
    # Append only the last character of next node
    reconstructed_string <- paste0(reconstructed_string, substr(next_node, nchar(next_node), nchar(next_node)))
  }
  
  # Remove first and last character since it's a cycle
  reconstructed_string <- substr(reconstructed_string, 2 , nchar(reconstructed_string) - (k-1))
  
  return(reconstructed_string)
}

k <- 8

kmers <- generate_binary_combinations(k)

# Generate overlap graph
result <- deBruijnGraph_test(kmers)

# Format the output as "KEY : VALUE1 VALUE2 ..."
result <- sapply(names(result), function(key) {
  paste0(key, ":", " ", paste(result[[key]], collapse = " "))
})

# Example edge order (output from Eulerian cycle function)
edge_order <- find_eulerian_cycle(result)  # Example cycle

# Reconstruct original string
reconstructed <- reconstruct_string(edge_order)
print(reconstructed)
