## nearest neighbour:

# Helper: null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

parse_edge_tree <- function(lines) {
  u <- lines[1]
  v <- lines[2]
  
  edges <- list()
  node_counts <- list()
  
  # Build edge list and count node occurrences
  for (i in 3:length(lines)) {
    parts <- strsplit(lines[i], "->")[[1]]
    parent <- parts[1]
    child  <- parts[2]
    edges[[length(edges) + 1]] <- c(parent, child)
    
    node_counts[[parent]] <- (node_counts[[parent]] %||% 0) + 1
    node_counts[[child]]  <- (node_counts[[child]] %||% 0) + 1
  }
  
  # Classify nodes for an unrooted tree
  classify_node <- function(node) {
    deg <- node_counts[[node]]
    if (deg == 2) {
      return("leaf")
    } else if (deg >= 6) {
      return("internal")
    } else {
      return("root") # degree 2 nodes shouldn't exist in a proper unrooted binary tree
    }
  }
  
  node_types <- sapply(names(node_counts), classify_node)
  
  list(edges = edges,u = u, v = v,
       leaves = names(node_types[node_types == "leaf"]),
       internal_nodes = names(node_types[node_types == "internal"]),
       root = names(node_types[node_types == "root"]))
}

dedup_edges_tree <- function(edges, leaves) {
  out <- list()
  seen <- character(0)
  
  is_leaf <- function(x) x %in% leaves
  
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

nearest_neighbors <- function(res) {
  u <- res$u
  v <- res$v
  internal <- res$internal_nodes
  leaves <- res$leaves
  edges <- res$edges
  
  if (!(u %in% internal && v %in% internal)) {
    stop("NNI requires u and v to be internal nodes")
  }
  
  # edges incident to u and v
  edges_u <- Filter(function(e) u %in% e, edges)
  edges_v <- Filter(function(e) v %in% e, edges)
  
  # central edge(s)
  edges_u_v <- Filter(function(e) (v %in% e), edges_u)
  
  # normalize
  edges_u <- dedup_edges_tree(edges_u, leaves)
  edges_v <- dedup_edges_tree(edges_v, leaves)
  
  # exclude central edge
  edges_u <- Filter(function(e) !(v %in% e), edges_u)
  edges_v <- Filter(function(e) !(u %in% e), edges_v)
  
  # other endpoints
  x_1 <- edges_u[[1]]; x_1 <- x_1[x_1 != as.character(u)]
  x_2 <- edges_u[[2]]; x_2 <- x_2[x_2 != as.character(u)]
  y_1 <- edges_v[[1]]; y_1 <- y_1[y_1 != as.character(v)]
  y_2 <- edges_v[[2]]; y_2 <- y_2[y_2 != as.character(v)]
  
  uc <- as.character(u)
  vc <- as.character(v)
  
  # Neighbor 1: swap (u,x_2) with (v,y_1)
  edges_u_1 <- lapply(edges_u, function(e) {
    if (uc %in% e && x_2 %in% e) {
      if (e[1] == uc && e[2] == x_2) c(uc, y_1) else c(y_1, uc)
    } else e
  })
  edges_v_1 <- lapply(edges_v, function(e) {
    if (vc %in% e && y_1 %in% e) {
      if (e[1] == vc && e[2] == y_1) c(vc, x_2) else c(x_2, vc)
    } else e
  })
  
  # Neighbor 2: swap (u,x_2) with (v,y_2)
  edges_u_2 <- lapply(edges_u, function(e) {
    if (uc %in% e && x_2 %in% e) {
      if (e[1] == uc && e[2] == x_2) c(uc, y_2) else c(y_2, uc)
    } else e
  })
  edges_v_2 <- lapply(edges_v, function(e) {
    if (vc %in% e && y_2 %in% e) {
      if (e[1] == vc && e[2] == y_2) c(vc, x_2) else c(x_2, vc)
    } else e
  })
  
  # Add reversed orientation
  add_reverse <- function(edges) {
    revs <- lapply(edges, function(e) c(e[2], e[1]))
    c(edges, revs)
  }
  
  edges_1 <- add_reverse(c(edges_u_1, edges_v_1))
  edges_2 <- add_reverse(c(edges_u_2, edges_v_2))
  
  # add central edge back
  edges_1 <- c(edges_1, edges_u_v)
  edges_2 <- c(edges_2, edges_u_v)
  
  # keep all other edges unchanged
  neighbour_cut <- Filter(function(e) !(u %in% e) && !(v %in% e), edges)
  
  neighbour_1 <- c(neighbour_cut, edges_1)
  neighbour_2 <- c(neighbour_cut, edges_2)
  
  list(neighbour1 = neighbour_1, neighbour2 = neighbour_2)
}

format_edges <- function(edges) {
  # edges is a list of character vectors c(u,v)
  lines <- vapply(edges, function(e) paste0(e[1], "->", e[2]), character(1))
  paste(lines, collapse = "\n")
}


# Example run
D <- "58 59
0->32
32->0
1->32
32->1
2->33
33->2
3->33
33->3
4->34
34->4
5->34
34->5
6->35
35->6
7->35
35->7
8->36
36->8
9->36
36->9
10->37
37->10
11->37
37->11
12->38
38->12
13->38
38->13
14->39
39->14
15->39
39->15
16->40
40->16
17->40
40->17
18->41
41->18
19->41
41->19
20->42
42->20
21->42
42->21
22->43
43->22
23->43
43->23
24->44
44->24
25->44
44->25
26->45
45->26
27->45
45->27
28->46
46->28
29->46
46->29
30->47
47->30
31->47
47->31
33->48
48->33
43->48
48->43
44->49
49->44
38->49
49->38
39->50
50->39
36->50
50->36
32->51
51->32
49->51
51->49
48->52
52->48
47->52
52->47
51->53
53->51
46->53
53->46
53->54
54->53
34->54
54->34
45->55
55->45
40->55
55->40
42->56
56->42
55->56
56->55
56->57
57->56
35->57
57->35
41->58
58->41
57->58
58->57
54->59
59->54
58->59
59->58
59->60
60->59
52->60
60->52
37->61
61->37
50->61
61->50
60->61
61->60"
lines <- unlist(strsplit(D," |\t|\n"))
res <- parse_edge_tree(lines)
neighbors <- nearest_neighbors(res)

cat(format_edges(neighbors$neighbour1), "\n\n", format_edges(neighbors$neighbour2),sep = "",file = "result.txt")
