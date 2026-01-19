### Amino trees
library(tidyverse)
input_file <- choose.files()
Amino_table <- read.csv(input_file,sep = " ",header = FALSE)
colnames(Amino_table) <- c("amino_letter","amino_trip","weight")

# Build AA_map: named list of lists
AA_map <- setNames(
  lapply(seq_len(nrow(Amino_table)), function(i) {
    list(
      Triplet = Amino_table$amino_trip[i],
      Mass    = Amino_table$weight[i]
    )
  }),
  Amino_table$amino_letter
)

spectrum_input <- "87 101 158 232 245 345 382 508 510 657 671 728 758 827 886 957 990 1028 1137 1159 1274 1290 1371 1387 1502 1524 1633 1671 1704 1775 1834 1903 1933 1990 2004 2151 2153 2279 2316 2416 2429 2503 2560 2574 2661"
spectrum_input <- sort(c(0,as.integer(unlist(strsplit(spectrum_input," ")))))

find_possible_path <- function(spectrum_input, Amino_table, AA_map) {
  # Step 1: define sink
  sink <- max(spectrum_input)
  
  # Step 2: construct graph edges
  edges <- character()
  for (i in seq_along(spectrum_input)) {
    for (j in seq((i+1), length(spectrum_input))) {
      diff <- spectrum_input[j] - spectrum_input[i]
      aa <- Amino_table[Amino_table[,"weight"]==diff,"amino_letter"]
      if (length(aa) > 0) {
        edges <- c(edges, sprintf("%d->%d:%s", spectrum_input[i], spectrum_input[j], aa[1]))
      }
    }
  }
  
  # Step 3: split into df
  parts <- strsplit(edges, "->|:")
  df <- do.call(rbind, lapply(parts, function(x) x))
  colnames(df) <- c("from", "to", "aa")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  Total_weight <- spectrum_input[[length(spectrum_input)]]
  
  # Step 4: recursive trace_back
  trace_back <- function(df, target) {
    matches <- df[df$to == target, ]
    if (nrow(matches) == 0) {
      return(list(list(nodes = target, aas = character())))
    }
    chains <- list()
    for (i in seq_len(nrow(matches))) {
      from_node <- matches$from[i]
      aa_label  <- matches$aa[i]
      preds <- trace_back(df, from_node)
      for (p in preds) {
        new_nodes <- c(p$nodes, target)
        new_aas   <- c(p$aas, aa_label)   # append at end
        chains <- c(chains, list(list(nodes = new_nodes, aas = new_aas)))
      }
    }
    chains
  }
  
  chains <- trace_back(df, Total_weight)
  
  # Step 5: collapse aas into strings and compute spectra
  rev_aas <- lapply(chains, function(ch) {
    ch$aas <- paste0(ch$aas, collapse = "")
    ch$spectra <- Peptide_linear_Spectrum(ch$aas, AA_map)
    ch
  })
  
  # Step 6: find peptide whose spectrum matches input
  answer_peptide <- character()
  for (i in seq_along(rev_aas)) {
    current_element <- rev_aas[[i]]
    current_spectrum <- current_element$spectra
    if (all(spectrum_input %in% current_spectrum)) {
      answer_peptide <- current_element$aas
      break
    }
  }
  
  return(answer_peptide)
}

answer <- find_possible_path(spectrum_input, Amino_table, AA_map)
print(answer)
