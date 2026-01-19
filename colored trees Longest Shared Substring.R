# colored trees Longest Shared Substring

# Parse adjacency list and leaf colors
parse_tree <- function(lines) {
  sep_index <- which(lines == "-")
  adj_lines <- lines[1:(sep_index-1)]
  color_lines <- lines[(sep_index+1):length(lines)]
  
  # adjacency list
  adj <- list()
  for (line in adj_lines) {
    parts <- strsplit(line, ":")[[1]]
    node <- trimws(parts[1])
    children <- character(0)
    if (length(parts) > 1 && nchar(trimws(parts[2])) > 0) {
      children <- strsplit(trimws(parts[2]), "\\s+")[[1]]
    }
    adj[[node]] <- children
  }
  
  # initial colors
  colors <- list()
  for (line in color_lines) {
    parts <- strsplit(line, "\\s+")[[1]]
    node <- parts[1]; col <- parts[2]
    colors[[node]] <- col
  }
  
  list(adj=adj, colors=colors)
}

# TreeColoring algorithm
tree_coloring <- function(adj, colors) {
  nodes <- names(adj)
  # initialize gray for unlabeled nodes
  for (n in nodes) {
    if (is.null(colors[[n]])) colors[[n]] <- "gray"
  }
  
  repeat {
    ripe_nodes <- c()
    for (n in nodes) {
      if (colors[[n]] == "gray") {
        children <- adj[[n]]
        if (length(children) > 0 && all(sapply(children, function(c) colors[[c]] != "gray"))) {
          ripe_nodes <- c(ripe_nodes, n)
        }
      }
    }
    if (length(ripe_nodes) == 0) break
    
    for (n in ripe_nodes) {
      child_colors <- unique(unlist(lapply(adj[[n]], function(c) colors[[c]])))
      if (length(child_colors) > 1) {
        colors[[n]] <- "purple"
      } else {
        colors[[n]] <- child_colors[1]
      }
    }
  }
  colors
}

# -------------------------------
# Example usage
# -------------------------------
# Convert a multi-line string into a character vector
to_lines <- function(block) {
  # Split on newline characters
  lines <- unlist(strsplit(block, "\n"))
  # Trim whitespace around each line
  lines <- trimws(lines)
  # Drop empty lines if any
  lines <- lines[lines != ""]
  return(lines)
}

# Example usage
block <- "0: 1 2
1: 3 4 5
2: 6 7 8
3: 25 26
4: 33 34
5: 9 10
6: 47 48 49
7: 16 17
8: 13 14 15
9:
10: 11 12
11:
12:
13:
14:
15:
16: 18 19 20
17: 27 28 29
18: 23 24
19: 21 22
20:
21:
22:
23: 30 31 32
24:
25:
26:
27: 35 36 37
28: 40 41 42
29:
30:
31: 38 39
32:
33:
34: 43 44
35:
36:
37:
38:
39:
40: 45 46
41:
42:
43:
44:
45:
46:
47:
48:
49:
-
9 blue
11 red
12 red
13 blue
14 blue
15 blue
20 blue
21 blue
22 red
24 blue
25 blue
26 blue
29 red
30 red
32 red
33 red
35 blue
36 red
37 red
38 red
39 red
41 red
42 red
43 red
44 red
45 blue
46 red
47 red
48 blue
49 red"

lines <- to_lines(block)
parsed <- parse_tree(lines)
result <- tree_coloring(parsed$adj, parsed$colors)

# Print in required format
for (n in sort(as.integer(names(result)))) {
  cat(n, result[[as.character(n)]], "\n")
}
