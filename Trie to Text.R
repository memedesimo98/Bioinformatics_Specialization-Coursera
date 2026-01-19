### Trie to Text

Trie_to_text <- function(text,trie_final) {
  
  splitted_text <- unlist(strsplit(text,""))
  Text_position <- 0
  max_position <- length(splitted_text) - 1
  while (Text_position <= max_position) {
    cat("\n--- Starting at Text_position:", Text_position, " ---\n")
    processing <- TRUE
    current_node_name <- "0"
    relative_text_position <- Text_position
    hard_fail <- FALSE
    while (processing && relative_text_position <= max_position) {
      relative_text_position <- relative_text_position + 1
      current_letter <- splitted_text[[relative_text_position]]
      current_node <- trie_final[[current_node_name]]
      current_node_name <- current_node$name
      
      cat("At node", current_node_name,
          "looking for letter", current_letter,
          "symbols:", current_node$symbol, "\n")
      
      if (!current_letter %in% current_node$symbol) {
        cat("Letter", current_letter, "not found at node", current_node_name, "\n")
        processing <- FALSE
        hard_fail <- TRUE
        break
      }
      
      pos <- which(current_node$symbol == current_letter)
      next_node <- current_node$children[[pos]]
      cat("Following edge", current_letter, "to node", next_node, "\n")
      
      if (!is.null(trie_final[[next_node]])) {
        current_node <- trie_final[[next_node]]
        current_node_name <- current_node$name
        cat("Moved to node", current_node_name, "\n")
      } else {
        cat("Next node", next_node, "does not exist\n")
        processing <- FALSE
        break
      }
    }
    if (!hard_fail) {
      if (!is.null(trie_final[[current_node_name]]$ending)) {
        cat("Pattern ended at node", current_node_name,
            "pattern:", trie_final[[current_node_name]]$ending, "\n")
        if (is.null(trie_final[[current_node_name]]$matches)) {
          trie_final[[current_node_name]]$matches <- Text_position
        } else {
          trie_final[[current_node_name]]$matches <- c(trie_final[[current_node_name]]$matches, Text_position)
        }
      }
    }
    
    Text_position <- Text_position + 1
  }
  return(trie_final)
}


text <- "AATCGGGTTCAATCGGGGT"

# Example: reading space-separated input and printing exactly as triples
input <- "CTTTCTGCT TATGCCTTA CCCGCGACC AAACGATAA TACATGTTA CTCTACCCT"
patterns <- strsplit(input, "\\s+", perl = TRUE)[[1]]

trie <- TrieConstruction(patterns)

test <- AnnotateTrieTerminals(trie,patterns)

final_test <- Trie_to_text(text,test)

for (i in 1:length(final_test)) {
  current_element <- final_test[[i]]
  if (!is.null(current_element$matches)) {
    cat(paste0(current_element$ending,":",collapse = ""), paste0(current_element$matches, collapse = " "), "\n")
  }
}

# test



