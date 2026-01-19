### RNA translation

# codon map

# Choose the file and store the path
file_path <- file.choose("RNA")

# Read the file
codon_data <- read.table(file_path, stringsAsFactors = FALSE)

# Choose the file and store the path
file_path <- file.choose("Amino")

# Read the file
AA_data <- read.table(file_path, stringsAsFactors = FALSE)

codon_data <- codon_data[codon_data[,2]!="-",]

AA_mapping <- function(codon_data,AA_data, codon_map) {
  
  # Create an empty list to store the mapping
  codon_map <- list()
  
  # Populate the mapping
  for (i in 1:nrow(codon_data)) {
    amino_acid <- codon_data[i, 2]  # Second column contains amino acid
    codon <- codon_data[i, 1]  # First column contains codon
    
    # If the amino acid doesn't exist in the list, initialize it as a named list
    if (!is.list(codon_map[[amino_acid]])) {
      codon_map[[amino_acid]] <- list(Codons = c())
    }
    
    # Append codon to the Codons list within the amino acid entry
    codon_map[[amino_acid]]$Codons <- c(codon_map[[amino_acid]]$Codons, codon)
  }
  
  Full_aa_map <- lapply(names(codon_map), function(aa) {
    # Extract the corresponding data from AA_data
    triplet_value <- AA_data[AA_data[,1] == aa, 2]
    mass_value <- AA_data[AA_data[,1] == aa, 3]
    
    # Add new fields to the existing structure
    codon_map[[aa]]$Triplet <- triplet_value
    codon_map[[aa]]$Mass <- mass_value
    
    return(codon_map[[aa]])  # Return the modified list element
  })
  
  names(Full_aa_map) <- names(codon_map)  # Preserve amino acid names
  
  return(Full_aa_map)
}

# Apply the function to update codon_map
AA_map_realistic <- AA_mapping(codon_data,AA_data, codon_map)

