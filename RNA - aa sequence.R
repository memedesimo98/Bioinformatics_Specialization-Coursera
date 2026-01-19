# RNA - aa sequence

split_codons <- function(sequence,start = 1,n = 3) {
  # Split the sequence into triplets
  finish<-start+2
  codons <- substring(sequence, seq(start, nchar(sequence), by = n), seq(finish, nchar(sequence), by = n))
  return(codons)
}

# Choose the file and store the path
file_path <- file.choose("input")

# Open and read the content of the text file
file_content <- readLines(file_path)

# Clean the text by removing unwanted characters (slashes, etc.)
cleaned_data <- gsub("\\\\", "", file_content)  # Remove backslashes
cleaned_data <- gsub("\"", "", cleaned_data)   # Remove extra quotes if any

dna_sequence_raw <- cleaned_data

# Example usage
rna_sequence <- "CCCCGUACGGAGAUGAAA"
codon_list <- split_codons(rna_sequence)

print(codon_list)

translate_codons <- function(codon_list, AA_map) {
  sapply(codon_list, function(codon) {
    # Find the amino acid corresponding to the codon
    aa <- names(AA_map)[sapply(AA_map, function(x) codon %in% x$Codons)]
    if (length(aa) == 0) return("?")  # Handle unknown codons
    return(aa)
  })
}

# Example usage
translated_sequence <- translate_codons(codon_list, AA_map)
print(names(translated_sequence))
unnamed_translated_sequence<-unname(translated_sequence)
print(paste(unnamed_translated_sequence, collapse = ""))


split_stop_codons <- function(named_vec, complete_list = TRUE) {
  # Find indices where NaN occurs
  nan_indices <- which(named_vec == "NaN")
  
  # Handle case when no NaN is found
  if (length(nan_indices) == 0) {
    return(list("1" = list(codon_sequence = names(named_vec), aa_sequence = as.character(named_vec))))
  }
  
  # Define break points based on complete_list argument
  break_points <- if (complete_list) c(0, nan_indices, length(named_vec) + 1) else c(0, nan_indices[1])
  
  # Apply function over the break points efficiently
  result_list <- setNames(
    lapply(seq_along(break_points[-1]), function(i) {
      start_idx <- break_points[i] + 1
      end_idx <- break_points[i + 1] - 1
      if (start_idx <= end_idx) {
        group <- named_vec[start_idx:end_idx]
        list(codon_sequence = names(group), aa_sequence = as.character(group))
      } else {
        NULL # Handles edge case where range is invalid
      }
    }), as.character(seq_along(break_points[-1])))
  
  # Remove NULL entries that might arise due to edge cases
  result_list <- result_list[!sapply(result_list, is.null)]
  
  return(result_list)
}


# Example usage
result <- split_stop_codons(translated_sequence)
print(result)

### works, next, from peptide to DNA
Find_sequence_from_peptide <- function(dna_sequence,peptide_sequence,AA_map) {
  
  result_list<-list()
  # convert reads into their complementary while mantaining their indexing
  convert_to_reverse_complementary <- function(dna_sequence) {
    base_complements <- c(A = "T", T = "A", C = "G", G = "C")
    parts <- strsplit(dna_sequence, "\\.")[[1]]
    main_sequence <- parts[1]
    index <- ifelse(length(parts) > 1, paste0(".", parts[2]), "")
    dna_bases <- unlist(strsplit(main_sequence, ""))
    complementary_bases <- base_complements[dna_bases]
    reverse_complementary_sequence <- paste(rev(complementary_bases), collapse = "")
    result <- paste0(reverse_complementary_sequence, index)
    return(result)
  }
  
  peptide_sequence<-unlist(strsplit(peptide_sequence,""))
  
  peptide_sequence_list<-lapply(peptide_sequence,function(x) AA_map[[x]]$Codons)
  
  DNA_to_RNA <- function(dna_sequence){
    dna_sequence<-unlist(strsplit(dna_sequence,""))
    dna_sequence[dna_sequence == "T"]<-"U"
    rna_sequence<-paste(dna_sequence,collapse="")
  }
  RNA_to_DNA <- function(rna_sequence){
    rna_sequence <- unlist(strsplit(rna_sequence, ""))
    rna_sequence[rna_sequence == "U"] <- "T"
    dna_sequence <- paste(rna_sequence, collapse = "")
    return(dna_sequence)
  }
  
  dna_complementary<-convert_to_reverse_complementary(dna_sequence)
  rna_sequence<-DNA_to_RNA(dna_sequence)
  rna_sequence_complementary<-DNA_to_RNA(dna_complementary)
  
  split_codons_frames <- function(sequence, N = 3) {
    frames <- list(
      "0" = split_codons(sequence, 1, N),
      "1" = split_codons(sequence, 2, N),
      "2" = split_codons(sequence, 3, N)
    )
    return(frames)
  }
  
  rna_codon_list<-split_codons_frames(rna_sequence)
  rna_complementary_codon_list<-split_codons_frames(rna_sequence_complementary)
  
  motif_codon_matches <- function(codon_sequence, peptide_sequence_list) {   
    codon_pos <- which(codon_sequence %in% peptide_sequence_list[[1]])  # Fix `which()` usage
    
    match_list <- list()
    
    if (length(codon_pos) > 0) {
      for (pos in codon_pos) {
        i <- 2
        starting_pos <- pos
        
        while (i <= length(peptide_sequence_list)) {
          match <- starting_pos + (i - 1)  # Next expected codon position
          if (match <= length(codon_sequence) && codon_sequence[match] %in% peptide_sequence_list[[i]]) {
            i <- i + 1
          } else {
            break  # Exit loop if the chain breaks
          }
        }
        
        # If all motifs matched
        if (i > length(peptide_sequence_list)) {
          start_match <- starting_pos * 3 - 2
          codons_sequence <- codon_sequence[starting_pos:(starting_pos + length(peptide_sequence_list) - 1)]
          match_list[[length(match_list) + 1]] <- list(Starting_pos = start_match,
                                                       codons = codons_sequence,
                                                       Sequence = RNA_to_DNA(paste(codons_sequence, collapse = "")))
        }
      }
    }
    
    return(match_list)
  }
  
  # Apply function across codon frames
  rna_match <- lapply(rna_codon_list, motif_codon_matches, peptide_sequence_list)
  rna_complementary_match <- lapply(rna_complementary_codon_list, motif_codon_matches, peptide_sequence_list)
  
  adjust_starting_pos <- function(named_list) {
    lapply(names(named_list), function(name) {
      increment <- as.numeric(name)  # Convert name to numeric for incrementing
      lapply(named_list[[name]], function(sub_list) {
        sub_list$Starting_pos <- sub_list$Starting_pos + increment
        return(sub_list)
      })
    })
  }
  
  rna_match<-adjust_starting_pos(rna_match)
  rna_complementary_match_raw<-adjust_starting_pos(rna_complementary_match)
  
  rna_complementary_match <- lapply(rna_complementary_match_raw, function(x) {
    lapply(x, function(y) {
      if (is.list(y) && "Sequence" %in% names(y)) { 
        y$Sequence <- convert_to_reverse_complementary(y$Sequence)
      }
      y
    })
  })
  
  
  
  result_list<-list(rna_match,rna_complementary_match)
  
  return(result_list)
}

dna_sequence<-"ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
peptide_sequence<-"MA"
result<-Find_sequence_from_peptide(dna_sequence,peptide_sequence,AA_map)

result_sequence <- unlist(lapply(result, function(x) {
  lapply(x, function(y) {
    lapply(y, function(z) {
      z$Sequence
    })
  })
}))
cat(paste(result_sequence,collapse = "\n"))

