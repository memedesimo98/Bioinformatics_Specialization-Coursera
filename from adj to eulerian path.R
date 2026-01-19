
# eulerian path function

# Function to parse adjacency list from formatted input
parse_adj_list <- function(input_lines) {
  adj_list <- list()
  
  for (line in input_lines) {
    parts <- unlist(strsplit(line, ":"))
    node <- trimws(parts[1])  # Extract node
    neighbors <- unlist(strsplit(trimws(parts[2]), " "))  # Extract neighbors
    
    adj_list[[node]] <- neighbors
  }
  
  return(adj_list)
}

# Function to find Eulerian path length & start node
find_eulerian_path_info <- function(adj_list) {
  out_degree <- list()
  in_degree <- list()
  
  # Calculate out-degree and in-degree
  for (node in names(adj_list)) {
    out_degree[[node]] <- length(adj_list[[node]])
    for (neighbor in adj_list[[node]]) {
      in_degree[[neighbor]] <- ifelse(is.null(in_degree[[neighbor]]), 1, in_degree[[neighbor]] + 1)
    }
  }
  
  # Compute total edges
  total_edges <- sum(unlist(out_degree))
  
  # Identify start node (out-degree = in-degree + 1)
  start_node <- NULL
  for (node in names(out_degree)) {
    in_deg <- ifelse(is.null(in_degree[[node]]), 0, in_degree[[node]])
    if (out_degree[[node]] == in_deg + 1) {
      start_node <- node
      break
    }
  }
  
  return(list("Path_Length" = total_edges + 1, "Start_Node" = start_node))
}

# Function to compute Eulerian path(s) using adj_list, Path_length, and Start_node
find_eulerian_paths <- function(adj_list, Path_length, Start_node) {
  # Initialize path exploration lists
  path_list <- list(list(
    current_node = Start_node,
    next_node = NA,
    path = c(Start_node),
    missing_length = Path_length - 1,  # Since path starts with Start_node
    node_adj_list = adj_list
  ))
  
  # Process paths iteratively
  complete_paths <- list()
  
  while (length(path_list) > 0) {
    new_paths <- list()
    
    for (entry in path_list) {
      # Extract the current state
      current_node <- entry$current_node
      node_adj_list <- entry$node_adj_list
      path <- entry$path
      missing_length <- entry$missing_length
      
      # Get possible neighbors for the current node
      if (!is.null(node_adj_list[[current_node]]) && length(node_adj_list[[current_node]]) > 0) {
        neighbors <- node_adj_list[[current_node]]
        
        for (next_node in neighbors) {
          # Create a new instance of the list for branching
          new_entry <- entry
          new_entry$next_node <- next_node
          new_entry$path <- c(path, next_node)
          new_entry$missing_length <- missing_length - 1
          
          # Remove used edge from adjacency list copy
          new_entry$node_adj_list[[current_node]] <- setdiff(neighbors, next_node)
          
          # Update current_node for next iteration
          new_entry$current_node <- next_node
          
          # If missing_length reaches zero, store as a valid path
          if (new_entry$missing_length == 0) {
            complete_paths[[length(complete_paths) + 1]] <- new_entry$path
          } else {
            new_paths[[length(new_paths) + 1]] <- new_entry
          }
        }
      }
    }
    
    # Update the path list with new branches
    path_list <- new_paths
  }
  
  return(complete_paths)
}

# Choose the file and store the path
file_path <- file.choose("text")

# Open and read the content of the text file
file_content <- readLines(file_path)

# Clean the text by removing unwanted characters (slashes, etc.)
cleaned_data <- gsub("\\\\", "", file_content)  # Remove backslashes
cleaned_data <- gsub("\"", "", cleaned_data)   # Remove extra quotes if any

kmers <- cleaned_data

data <- data[order(sapply(data, function(line) as.numeric(strsplit(line, ":")[[1]][1])))]

kmers<-unlist(strsplit(kmers," "))

# Function to generate all binary strings of length k
generate_binary_combinations <- function(k) {
  # Create all possible combinations using expand.grid()
  binary_matrix <- expand.grid(rep(list(c("0", "1")), k))
  
  # Convert matrix rows into character strings
  binary_strings <- apply(binary_matrix, 1, paste0, collapse = "")
  
  return(binary_strings)
}

k <- 3

kmers <- generate_binary_combinations(k)

kmers<-convert_binary_to_letters(kmers)

adj_list <- deBruijnGraph_test(kmers)

adj_list <- parse_adj_list(kmers)
# Compute Path_length and Start_node using previous function
eulerian_info <- find_eulerian_path_info(adj_list)

eulerian_info$Start_Node<-"AAA"
# Execute Eulerian path computation
paths <- find_eulerian_paths(adj_list, eulerian_info$Path_Length, eulerian_info$Start_Node)

# Output found paths
print(paste(paths[[1]], collapse = " "))

# Function to extract the last letter of a given string
last_letter <- function(word) {
  return(substr(word, nchar(word), nchar(word)))
}
first_entry<-paths[[1]][1]
rest_of_list<-paths[[1]][-1]
final_path<-paste0(first_entry,paste(sapply(rest_of_list,last_letter), collapse = ""),collapse = "")
final_path



# Function to transform 0 -> A and 1 -> B
convert_binary_to_letters <- function(binary_vector) {
  return(chartr("01", "AB", binary_vector))
}

# Function to revert back from A -> 0 and B -> 1
convert_letters_to_binary <- function(letter_vector) {
  return(chartr("AB", "01", letter_vector))
}

# Example usage
binary_input <- c("000", "101", "111", "010")
converted_letters <- convert_binary_to_letters(binary_input)
original_binary <- convert_letters_to_binary(converted_letters)

print(converted_letters)  # Output: c("AAA", "BAB", "BBB", "ABA")
print(original_binary)    # Output: c("000", "101", "111", "010")
