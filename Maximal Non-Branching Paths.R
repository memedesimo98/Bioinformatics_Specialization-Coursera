### Maximal Non-Branching Paths 

parse_input <- function(lines) {
  edges <- list()
  leaves <- list()
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[i], ":")[[1]]
    parent <- trimws(parts[1])
    if (length(parts) > 1) {
      children <- strsplit(trimws(parts[2]), "\\s+")[[1]]
      for (child in children) {
        edges[[length(edges) + 1]] <- c(parent, child)
        if (grepl("^[ACGT]+$", child)) {
          leaves[[child]] <- child
        }
      }
    }
  }
  list(edges = edges, leaves = leaves)
}


children_map <- function(edges) {
  ch <- list()
  for (e in edges) {
    p <- e[1]; c <- e[2]
    ch[[p]] <- c(ch[[p]], c)
  }
  ch
}

all_nodes <- function(edges) {
  unique(unlist(edges))
}

# Degrees
compute_degrees <- function(edges) {
  nodes <- all_nodes(edges)
  indeg <- setNames(rep(0, length(nodes)), nodes)
  outdeg <- setNames(rep(0, length(nodes)), nodes)
  
  for (e in edges) {
    p <- e[1]; c <- e[2]
    outdeg[p] <- outdeg[p] + 1
    indeg[c] <- indeg[c] + 1
  }
  list(indeg=indeg, outdeg=outdeg)
}

# Maximal NonBranching Paths
maximal_nonbranching_paths <- function(edges) {
  ch <- children_map(edges)
  degs <- compute_degrees(edges)
  paths <- list()
  visited <- list()
  
  # Paths starting from nodes not 1-in-1-out
  for (node in names(ch)) {
    indeg <- degs$indeg[node]; outdeg <- degs$outdeg[node]
    if (!(indeg == 1 && outdeg == 1)) {
      if (outdeg > 0) {
        for (nbr in ch[[node]]) {
          path <- c(node, nbr)
          current <- nbr
          while (degs$indeg[current] == 1 && degs$outdeg[current] == 1) {
            nextnode <- ch[[current]][1]
            path <- c(path, nextnode)
            current <- nextnode
          }
          paths <- append(paths, list(path))
        }
      }
    }
  }
  
  # Cycles
  for (node in names(ch)) {
    indeg <- degs$indeg[node]
    outdeg <- degs$outdeg[node]
    if (indeg == 1 && outdeg == 1 && is.null(visited[[node]])) {
      # Only proceed if this node has exactly one child
      if (length(ch[[node]]) == 1) {
        cycle <- c(node)
        current <- ch[[node]][1]
        while (current != node) {
          # stop if current doesn't have exactly one child
          if (length(ch[[current]]) != 1) break
          cycle <- c(cycle, current)
          visited[[current]] <- TRUE
          current <- ch[[current]][1]
        }
        # only accept if we really closed the cycle
        if (current == node) {
          cycle <- c(cycle, node)
          paths <- append(paths, list(cycle))
        }
      }
    }
  }
  
  
  paths
}

lines <- c("1: 2", "2: 3", "3: 4 5", "6: 7", "7: 6")
parsed <- parse_input(lines)
paths <- maximal_nonbranching_paths(parsed$edges)

for (p in paths) {
  cat(paste(p, collapse=" "), "\n")
}

