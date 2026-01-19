
parse_input <- function(lines) {
  n <- as.integer(lines[1])
  edges <- list()
  leaves <- list()
  for (i in 2:length(lines)) {
    parts <- strsplit(lines[i], "->")[[1]]
    parent <- parts[1]
    child  <- parts[2]
    edges[[length(edges) + 1]] <- c(parent, child)
    if (grepl("^[ACGT]+$", child)) {
      # Leaf node is the DNA string itself
      leaves[[child]] <- child
    }
  }
  list(n = n, edges = edges, leaves = leaves)
}

hamming <- function(a, b) {
  sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
}

adjacency_from_chmap <- function(chmap, root, leaves, internal_nodes,bidirectional = TRUE) {
  if (bidirectional) {
    bid_adj <- list()
  }
  adj <- list()
  for (p in names(chmap)) {
    if (p == root) next
    label_p <- if (p %in% names(internal_nodes)) internal_nodes[[p]] else leaves[[p]]
    for (u in chmap[[p]]) {
      if (u == root) next
      label_u <- if (u %in% names(internal_nodes)) internal_nodes[[u]] else leaves[[u]]
      d <- hamming(label_p, label_u)
      adj[[length(adj) + 1]] <- paste0(label_p, "->", label_u, ":", d)
      if (bidirectional) {
        bid_adj[[length(bid_adj) + 1]] <- paste0(label_p, "->", label_u, ":", d)
        bid_adj[[length(bid_adj) + 1]] <- paste0(label_u, "->", label_p, ":", d)
      }
    }
  }
  if (bidirectional) {
    tot_adj <- list(adj,bid_adj)
    return(tot_adj)
  } else {
    adj
    }
}



# Build a directed children map starting from a root node
build_children_map <- function(edges, root, leaves) {
  # Build adjacency (undirected)
  adj <- lapply(unique(unlist(edges)), function(v) character(0))
  names(adj) <- unique(unlist(edges))
  for (e in edges) {
    adj[[e[1]]] <- c(adj[[e[1]]], e[2])
    adj[[e[2]]] <- c(adj[[e[2]]], e[1])
  }
  
  # DFS/BFS to orient edges away from root
  chmap <- list()
  visited <- c(root)
  queue <- c(root)
  while (length(queue) > 0) {
    v <- queue[1]; queue <- queue[-1]
    kids <- setdiff(adj[[v]], visited)
    if (length(kids) > 0) {
      chmap[[v]] <- kids
      visited <- c(visited, kids)
      queue <- c(queue, kids)
    }
  }
  chmap
}


# Build a parallel chmap: leaves become singleton sets, internals stay as IDs
init_parallel_chmap <- function(chmap, leaf_chars) {
  is_internal <- function(x) x %in% names(chmap)
  out <- lapply(chmap, function(kids) {
    lapply(kids, function(u) if (is_internal(u)) u else leaf_chars[[u]])
  })
  names(out) <- names(chmap)
  out
}

# Replace every occurrence of an internal node ID with its computed set in the parallel chmap
propagate_set_to_parallel <- function(chmap_par, node_id, set_vec) {
  for (p in names(chmap_par)) {
    for (i in seq_along(chmap_par[[p]])) {
      child <- chmap_par[[p]][[i]]
      # If child is exactly the node ID string, replace with the set vector
      if (is.character(child) && length(child) == 1 && child == node_id) {
        chmap_par[[p]][[i]] <- set_vec
      }
    }
  }
  chmap_par
}

# Children are ready if both are letter sets (character vectors of alphabet symbols), not node IDs
children_ready <- function(kids_par) {
  all(vapply(kids_par, function(x) is.character(x) && length(x) >= 1 && !(length(x) == 1 && x %in% names(kids_par)), logical(1)))
}

# Fitch bottom-up over parallel chmap
fitch_bottom_up_queue <- function(chmap, leaf_chars) {
  # Initialize parallel map with leaves set
  chmap_par <- init_parallel_chmap(chmap, leaf_chars)
  
  sets  <- list()
  score <- 0
  queue <- names(chmap)  # internal nodes only
  
  # Process until all internals get sets
  while (length(queue) > 0) {
    v <- queue[1]; queue <- queue[-1]
    if (!is.null(sets[[v]])) next  # already processed
    
    kids_par <- chmap_par[[v]]     # two elements; each either a set (vector) or an ID string
    ready <- all(vapply(kids_par, function(x) is.character(x) && !(length(x) == 1 && x %in% names(chmap)), logical(1)))
    
    if (ready) {
      s1 <- kids_par[[1]]
      s2 <- kids_par[[2]]
      inter <- intersect(s1, s2)
      if (length(inter) > 0) {
        sets[[v]] <- inter
      } else {
        sets[[v]] <- sort(union(s1, s2))
        score <- score + 1
      }
      # propagate v's set upward so parents see ready children next time
      chmap_par <- propagate_set_to_parallel(chmap_par, v, sets[[v]])
    } else {
      queue <- c(queue, v)  # not ready yet; revisit later
    }
  }
  
  list(sets = sets, score = score)
}

fitch_top_down_assign <- function(root, chmap, sets) {
  assign <- list()
  dfs <- function(node, parent_state = NULL) {
    if (!(node %in% names(chmap))) return() # leaf, skip
    if (is.null(parent_state)) {
      state <- sort(sets[[node]])[1]
    } else if (parent_state %in% sets[[node]]) {
      state <- parent_state
    } else {
      state <- sort(sets[[node]])[1]
    }
    assign[[node]] <<- state
    for (u in chmap[[node]]) {
      if (u %in% names(chmap)) dfs(u, state)
    }
  }
  dfs(root)
  assign
}

# Function to separate internal edge and split graph
split_internal_edge <- function(edges, leaves) {
  # Convert leaves to a simple character vector
  leaf_nodes <- unlist(leaves, use.names = FALSE)
  
  # Find candidate edges (both ends not in leaves)
  internal_edges <- Filter(function(e) !(e[1] %in% leaf_nodes || e[2] %in% leaf_nodes), edges)
  
  if (length(internal_edges) == 0) {
    stop("No internal edges found")
  }
  
  # Take the first internal edge (you can generalize later)
  chosen <- internal_edges[[1]]
  node1 <- chosen[1]
  node2 <- chosen[2]
  
  # Remove chosen edge from edges
  remaining_edges <- Filter(function(e) !(all(e == chosen)), edges)
  
  # Build adjacency map
  adj <- lapply(unique(unlist(remaining_edges)), function(v) character(0))
  names(adj) <- unique(unlist(remaining_edges))
  for (e in remaining_edges) {
    adj[[e[1]]] <- c(adj[[e[1]]], e[2])
    adj[[e[2]]] <- c(adj[[e[2]]], e[1])
  }
  
  # BFS to collect connected components
  bfs <- function(start) {
    visited <- c()
    queue <- c(start)
    while (length(queue) > 0) {
      v <- queue[1]; queue <- queue[-1]
      if (!(v %in% visited)) {
        visited <- c(visited, v)
        queue <- c(queue, adj[[v]])
      }
    }
    visited
  }
  
  comp1 <- bfs(node1)
  comp2 <- bfs(node2)
  
  # Partition edges
  edges1 <- Filter(function(e) all(e %in% comp1), remaining_edges)
  edges2 <- Filter(function(e) all(e %in% comp2), remaining_edges)
  
  list(
    removed_edge = chosen,
    component1 = edges1,
    component2 = edges2
  )
}


dedup_edges <- function(edges) {
  out <- list()
  seen <- character(0)
  
  is_leaf <- function(x) grepl("^[ACGT-]+$", x)
  
  for (e in edges) {
    u <- e[1]; v <- e[2]
    key <- paste(sort(c(u,v)), collapse = "-")
    
    if (key %in% seen) next
    seen <- c(seen, key)
    
    # Rule 1: if one endpoint is a leaf, keep orientation parent->leaf
    if (is_leaf(v) && !is_leaf(u)) {
      out <- c(out, list(c(u,v)))
    } else if (is_leaf(u) && !is_leaf(v)) {
      out <- c(out, list(c(v,u)))
    } else {
      # Rule 2: both internal, keep orientation with smaller numeric as child
      nu <- suppressWarnings(as.numeric(u))
      nv <- suppressWarnings(as.numeric(v))
      if (!is.na(nu) && !is.na(nv)) {
        if (nu < nv) {
          out <- c(out, list(c(v,u)))  # smaller as child
        } else {
          out <- c(out, list(c(u,v)))
        }
      } else {
        # fallback: just keep u->v
        out <- c(out, list(c(u,v)))
      }
    }
  }
  out
}

small_parsimony_unrooted <- function(lines,bidirectional = TRUE) {
  parsed <- parse_input(lines)
  edges  <- dedup_edges(parsed$edges)
  leaves <- parsed$leaves
  alphabet <- c("A","C","G","T")
  
  # Split an internal edge and orient subtrees outward
  total_tree <- split_internal_edge(edges, leaves)
  subtree_u <- total_tree$component1
  subtree_v <- total_tree$component2
  
  chmap_u <- build_children_map(subtree_u, total_tree$removed_edge[1], leaves)
  chmap_v <- build_children_map(subtree_v, total_tree$removed_edge[2], leaves)
  
  # Create artificial root and reunite under a single directed map
  root <- as.character(max(as.numeric(names(chmap_u)), 
                           as.numeric(names(chmap_v))) + 1)
  chmap <- c(chmap_u, chmap_v)
  chmap[[root]] <- c(total_tree$removed_edge[1], total_tree$removed_edge[2])
  
  old_edge <- total_tree$removed_edge
  
  # Internal nodes are keys of chmap; leaves are the DNA strings appearing as children
  parents <- names(chmap)
  
  # String length
  m <- nchar(leaves[[1]])
  
  # Storage for reconstructed labels (internal only)
  internal_nodes <- setNames(rep("", length(parents)), parents)
  
  total_score <- 0
  
  for (pos in 1:m) {
    # Column characters for leaves (keyed by leaf DNA strings)
    leaf_chars <- lapply(leaves, function(s) substr(s, pos, pos))
    
    # Fitch bottom-up via queue over chmap
    fb <- fitch_bottom_up_queue(chmap, leaf_chars)
    sets <- fb$sets
    total_score <- total_score + fb$score
    
    # Fitch top-down assignment over chmap
    assign <- fitch_top_down_assign(root, chmap, sets)
    
    # Append assigned symbol to internal node labels only
    for (v in parents) {
      internal_nodes[[v]] <- paste0(internal_nodes[[v]], assign[[v]])
    }
  }
  
  # Your internal_root handling (kept as you want it)
  internal_root <- internal_nodes[setdiff(names(internal_nodes), root)]
  internal_root <- internal_nodes[!internal_nodes %in% internal_root]
  
  # Build adjacency from chmap (skip artificial root)
  total_adjacency <- adjacency_from_chmap(chmap, root, leaves, internal_nodes,bidirectional)
  if (bidirectional) {
    adjacency <- total_adjacency[[1]]
    bidirectional_adjacency <- total_adjacency[[2]]
  } else {
    adjacency <- total_adjacency
  }
  # Reinsert original broken edge with both directions
  u <- old_edge[1]; v <- old_edge[2]
  label_u <- internal_nodes[[u]]
  label_v <- internal_nodes[[v]]
  d <- hamming(label_u, label_v)
  adjacency[[length(adjacency) + 1]] <- paste0(label_u, "->", label_v, ":", d)
  if (bidirectional) {
    bidirectional_adjacency[[length(bidirectional_adjacency) + 1]] <- paste0(label_u, "->", label_v, ":", d)
    bidirectional_adjacency[[length(bidirectional_adjacency) + 1]] <- paste0(label_v, "->", label_u, ":", d)
  }
  if (bidirectional) {
    list(score = total_score, adjacency = adjacency,
         bidirectional_adjacency = bidirectional_adjacency,
         internal = internal_nodes, leaves = leaves, lines = lines)
  } else {
  list(score = total_score, adjacency = adjacency, internal = internal_nodes, leaves = leaves, lines = lines)
  }
}

res <- small_parsimony_unrooted(lines)
cat(res$score, paste(res$adjacency, collapse = "\n"), sep = "\n",
    file = "result.txt")
