# Define the unique single-character identifiers

stuff <- "ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõöøùúûüýþÿĀāĂăĄąĆćĈĉĊċČčĎďĐđĒēĔĕĖėĘęĚěĜĝĞğĠġĢģĤĥĦħĨĩĪīĬĭĮįİıĲĳĴĵĶķĸĹĺĻļĽľĿŀŁłŃńŅņŇňŉŊŋŌōŎŏŐőŒœŔŕŖŗŘřŚśŜŝŞşŠšŢţŤťŦŧŨũŪūŬŭŮůŰűŲųŴŵŶŷŸŹźŻżŽž"

keys <- unique(unlist(strsplit(stuff, "")))

keys <- keys[1:144]

# Assign Mass values from 57 to 200
masses <- as.numeric(seq(57, 200))

# Create named list of lists
mass_list <- setNames(lapply(masses, function(m) list(Mass = m)), keys)

# Print the list
print(mass_list)

str(AA_map)
str(mass_list)
N <- 1000

spectrum <- c(0, 97, 99, 114, 128, 147, 147, 163, 186, 227, 241, 242, 244, 260, 261, 262, 283, 291,
              333, 340, 357, 385, 389, 390, 390, 405, 430, 430, 447, 485, 487, 503, 504, 518, 543,
              544, 552, 575, 577, 584, 632, 650, 651, 671, 672, 690, 691, 738, 745, 747, 770, 778,
              779, 804, 818, 819, 820, 835, 837, 875, 892, 917, 932, 932, 933, 934, 965, 982, 989,
              1030, 1039, 1060, 1061, 1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194,
              1194, 1208, 1209, 1223, 1225, 1322)

result <- subpeptides_scoring_finder_ultimate(spectrum,N,mass_list)
result_final <- subpeptide_to_masses_reducted(result,mass_list)
paste(unique(result_final),collapse = " ")
