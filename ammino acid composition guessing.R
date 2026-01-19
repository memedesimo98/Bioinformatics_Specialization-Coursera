### ammino acid composition guessing

mass_guessing <- function(spectrum) {
sorted_spectrum <- rev(sort(spectrum))
maximum_length <- length(sorted_spectrum)

ghost_masses <- numeric()
for (mass in 1:maximum_length) {
  current_mass <- sorted_spectrum[mass]
  for (iteration in mass:(maximum_length)) {
    ghost_mass <- as.character(current_mass - sorted_spectrum[iteration])
    if (ghost_mass > 0) {
    ghost_masses <- c(ghost_masses,ghost_mass)
    }
  }
}
table_ghost_masses <- as.data.frame(table(ghost_masses))
sorted_table <- table_ghost_masses[order(-table_ghost_masses[,2]),]
sorted_ghost_masses <- character()
sorted_ghost_repetitions <- character()
for (i in 1:nrow(sorted_table)) {
  sorted_ghost_masses <- c(sorted_ghost_masses, rep(as.character(sorted_table[i,1]),sorted_table[i,2]))
  sorted_ghost_repetitions <- c(sorted_ghost_repetitions, rep(as.character(sorted_table[i,2]),sorted_table[i,2]))
}
sorted_ghost <- setNames(sorted_ghost_masses,sorted_ghost_repetitions)
return(sorted_ghost)
}

filter_aa_map <- function(mass_guess,M,AA_map) {
  names_vec <- unique(as.numeric(names(mass_guess)))
  if (M <= length(names_vec)) {
    names_vec <- names_vec[M]
  } else {
    names_vec <- names_vec[length(names_vec)]
  }
  filtered_mass_guess <- unique(mass_guess[as.numeric(names(mass_guess)) >= names_vec])
  filtered_mass_guess <- filtered_mass_guess[
    as.numeric(filtered_mass_guess) >= 57 & as.numeric(filtered_mass_guess) <= 200
  ]
  AA_map_filtered <- Filter(function(x) x$Mass %in% as.numeric(filtered_mass_guess), AA_map)
  return(AA_map_filtered)
}

spectrum <- c(0,57,118,179,236,240,301 )
M <- 20
N <- 1000
result <- mass_guessing(spectrum)
filtered_mass_list <- filter_aa_map(result,M,mass_list)
result_masses <- subpeptides_scoring_finder_ultimate(spectrum,N,filtered_mass_list)
result_final <- subpeptide_to_masses_reducted(result_masses,filtered_mass_list)
paste(result_final,collapse = " ")
