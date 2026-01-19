### Trie Construction

TrieConstruction <- function(patterns) {
  
  # Root node
  next_id <- 1L
  trie_final <- list("0" = list(name = "0",children = c(),symbol = c()))
  
  for (pattern in patterns) {
    current_node <- trie_final[["0"]]
    current_pattern <- unlist(strsplit(pattern,""))
    while (length(current_pattern)>0) {
      letter <- current_pattern[[1]]
      current_pattern <- current_pattern[-1]
      if (letter %in% current_node$symbol) {
        pos <- which(current_node$symbol == letter)
        next_node <- current_node$children[[pos]]
        current_node <- trie_final[[next_node]]
      } else {
        initial_next_id <- next_id
        current_node_name <- current_node$name
        while (length(current_pattern)>0) {
          trie_final[[current_node_name]]$children <- c(trie_final[[current_node_name]]$children,
                                                        as.character(next_id))
          trie_final[[current_node_name]]$symbol <- c(trie_final[[current_node_name]]$symbol,letter)
          letter <- current_pattern[[1]]
          current_pattern <- current_pattern[-1]
          current_node_name <- as.character(next_id)
          trie_final[[current_node_name]] <- list(name = current_node_name,children = c(),symbol = c())
          next_id <- next_id + 1
        }
        trie_final[[current_node_name]]$children <- c(trie_final[[current_node_name]]$children,
                                                      as.character(next_id))
        trie_final[[current_node_name]]$symbol <- c(trie_final[[current_node_name]]$symbol,letter)
        next_id <- next_id + 1
      }
    }
  }
  return(trie_final)
}

AnnotateTrieTerminals <- function(trie_final, patterns) {
  for (pattern in patterns) {
    current_node <- trie_final[["0"]]
    current_pattern <- unlist(strsplit(pattern,""))
    current_node_name <- current_node$name
    node_path <- c()
    for (i in 1:(length(current_pattern)-1)) {
      node_path <- c(node_path,current_node_name)
      current_letter <- current_pattern[[i]]
      pos <- which(current_node$symbol == current_letter)
      next_node <- current_node$children[[pos]]
      current_node <- trie_final[[next_node]]
      current_node_name <- current_node$name
    }
    node_path <- c(node_path,current_node_name)
    trie_final[[current_node_name]]$ending <- pattern
    trie_final[[current_node_name]]$node_path <- node_path
  }
  return(trie_final)
}


format_trie <- function(trie_final) {
  # iterate over all nodes in the trie
  for (node in trie_final) {
    # for each child edge, print "from to symbol"
    for (i in seq_along(node$children)) {
      cat(node$name, node$children[[i]], node$symbol[[i]], "\n")
    }
  }
}


# Example: reading space-separated input and printing exactly as triples
input <- "ATAGA ATC GAT"
patterns <- strsplit(input, "\\s+", perl = TRUE)[[1]]

adj <- TrieConstruction(patterns)

# Print adjacency list
format_trie(adj)

test <- AnnotateTrieTerminals(adj,patterns)

input <- choose.files()
input <- readLines(input)

# Write adjacency list to a text file outside the function
con <- file("trie_output.txt", open = "wt")

for (node in adj) {
  for (i in seq_along(node$children)) {
    line <- paste(node$name, node$children[[i]], node$symbol[[i]])
    writeLines(line, con)
  }
}

close(con)
