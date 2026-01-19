## Large parsimony, neighbour + unrooted

Large_unrooted_parsimony <- function(lines) {
  
  res <- small_parsimony_unrooted(lines)
  
  old_score <- res$score
  best_score <- old_score
  old_adj <- res$bidirectional_adjacency
  best_bidi_adj <- old_adj
  best_lines <- res$lines
  # defy n
  n <- length(res$leaves)
  process_finished <- FALSE
  
  while (process_finished == FALSE) {
    New_best_found <- FALSE
    
  # Collapse into one string with \n separators
  new_lines <- best_lines[-1]
  
  cleaned_edges <- strsplit(new_lines,"->")
  
  # Vector of leaf labels
  leaf_labels <- names(res$leaves)
  
  # Keep only edges where neither endpoint is a leaf
  filtered_edges <- Filter(function(e) {
    !(e[1] %in% leaf_labels || e[2] %in% leaf_labels)
  }, cleaned_edges)
  
  for (j in seq_along(filtered_edges)) {
    current_edge <- filtered_edges[[j]]
  new_lines_edged <- c(current_edge[[1]],current_edge[[2]],new_lines)
  resulting_trees <- parse_edge_tree(new_lines_edged)
  neighbours_tree <- nearest_neighbors(resulting_trees)
  neighbours_tree_1 <- format_edges(neighbours_tree$neighbour1)
  neighbours_tree_1 <- paste(n,neighbours_tree_1,sep = "\n")
  neighbours_tree_2 <- format_edges(neighbours_tree$neighbour2)
  neighbours_tree_2 <- paste(n,neighbours_tree_2,sep = "\n")
  neighbours_tree_lines <- list(neighbours_tree_1,neighbours_tree_2)
  for (i in 1:2) {
    temporary_lines <- neighbours_tree_lines[[i]]
    temporary_lines <- unlist(strsplit(temporary_lines," |\t|\n"))
    temporary_res <- small_parsimony_unrooted(temporary_lines)
    temporary_score <- temporary_res$score
    if (temporary_score<best_score) {
      old_score <- best_score
      old_adj <- best_bidi_adj
      best_score <- temporary_score
      best_bidi_adj <- temporary_res$bidirectional_adjacency
      best_lines <- temporary_res$lines
      New_best_found <- TRUE
    }
  }
  if (New_best_found) {
    break
  } else {
    if (j == length(filtered_edges)) {
      process_finished <- TRUE
  }
}
  }
  }
  return(list(old_score = old_score, old_adj = old_adj, best_score = best_score, best_bidi_adj = best_bidi_adj))
}


### output test with all best trees

Large_unrooted_parsimony_all_output <- function(lines) {
  
  res <- small_parsimony_unrooted(lines)
  
  old_score <- res$score
  best_score <- old_score
  old_adj <- res$bidirectional_adjacency
  best_bidi_adj <- old_adj
  best_lines <- res$lines
  n <- length(res$leaves)
  process_finished <- FALSE
  
  # container to collect outputs
  steps <- list()
  
  while (!process_finished) {
    New_best_found <- FALSE
    
    new_lines <- best_lines[-1]
    cleaned_edges <- lapply(new_lines, function(x) strsplit(x, "->")[[1]])
    leaf_labels <- names(res$leaves)
    
    filtered_edges <- Filter(function(e) {
      !(e[1] %in% leaf_labels || e[2] %in% leaf_labels)
    }, cleaned_edges)
    
    for (j in seq_along(filtered_edges)) {
      current_edge <- filtered_edges[[j]]
      new_lines_edged <- c(current_edge[[1]], current_edge[[2]], new_lines)
      resulting_trees <- parse_edge_tree(new_lines_edged)
      neighbours_tree <- nearest_neighbors(resulting_trees)
      
      neighbours_tree_lines <- list(
        paste(n, format_edges(neighbours_tree$neighbour1), sep="\n"),
        paste(n, format_edges(neighbours_tree$neighbour2), sep="\n")
      )
      
      for (i in 1:2) {
        temporary_lines <- neighbours_tree_lines[[i]]
        temporary_lines <- unlist(strsplit(temporary_lines," |\t|\n"))
        temporary_lines <- temporary_lines[nzchar(temporary_lines)]
        
        temporary_res <- small_parsimony_unrooted(temporary_lines)
        temporary_score <- temporary_res$score
        
        
        if (temporary_score < best_score) {
          old_score <- best_score
          old_adj <- best_bidi_adj
          best_score <- temporary_score
          best_bidi_adj <- temporary_res$bidirectional_adjacency
          best_lines <- temporary_res$lines
          New_best_found <- TRUE
          # record this neighbor step
          steps[[length(steps)+1]] <- paste(
            temporary_score,
            paste(
              # remove any line equal to n
              temporary_res$bidirectional_adjacency[temporary_res$bidirectional_adjacency != as.character(n)],
              collapse = "\n"
            ),
            sep = "\n"
          )
        }
      }
      
      if (New_best_found) {
        break
      } else if (j == length(filtered_edges)) {
        process_finished <- TRUE
      }
    }
  }
  
  # join steps with blank lines
  output <- paste(unlist(steps), collapse="\n\n")
  cat(output)
  
  invisible(list(
    old_score = old_score,
    old_adj = old_adj,
    best_score = best_score,
    best_bidi_adj = best_bidi_adj,
    steps = steps
  ))
}

result <- Large_unrooted_parsimony_all_output(lines)
output_text <- paste(unlist(result$steps), collapse = "\n\n")
writeLines(output_text, "NNI_output.txt")

