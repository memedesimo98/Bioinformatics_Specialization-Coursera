install.packages("igraph")

library(igraph)

# Sample input
input <- "32
60->61:6
60->11:6
60->59:6
61->60:6
61->18:5
61->57:15
18->61:5
57->61:15
57->14:14
57->56:6
11->60:6
59->60:6
59->12:6
59->58:6
12->59:6
58->59:6
58->8:5
58->50:12
8->58:5
50->58:12
50->1:6
50->49:14
14->57:14
56->57:6
56->6:9
56->55:6
6->56:9
55->56:6
55->24:12
55->54:11
24->55:12
54->55:11
54->30:10
54->53:7
30->54:10
53->54:7
53->4:13
53->52:10
4->53:13
52->53:10
52->22:13
52->51:7
22->52:13
51->52:7
51->15:14
51->43:13
15->51:14
43->51:13
43->19:12
43->41:13
1->50:6
49->50:14
49->17:11
49->48:11
17->49:11
48->49:11
48->26:8
48->47:6
26->48:8
47->48:6
47->21:14
47->46:6
21->47:14
46->47:6
46->10:15
46->45:9
10->46:15
45->46:9
45->9:7
45->44:11
9->45:7
44->45:11
44->28:9
44->42:12
28->44:9
42->44:12
42->3:8
42->39:13
19->43:12
41->43:13
41->5:10
41->40:11
3->42:8
39->42:13
39->27:11
39->38:10
5->41:10
40->41:11
40->16:9
40->20:13
16->40:9
20->40:13
27->39:11
38->39:10
38->7:10
38->37:15
7->38:10
37->38:15
37->25:9
37->36:6
25->37:9
36->37:6
36->23:6
36->35:10
23->36:6
35->36:10
35->29:11
35->34:9
29->35:11
34->35:9
34->2:9
34->33:12
2->34:9
33->34:12
33->13:8
33->32:7
13->33:8
32->33:7
32->0:11
32->31:12
0->32:11
31->32:12"

n <- 9
input <- unlist(strsplit(input,"\n"))
input <- formatted_new
# Parse number of leaves
n_leaves <- n
edges_raw <- input[-1]

# Parse edges into a data frame
parse_edge <- function(line) {
  parts <- strsplit(line, "->|:")[[1]]
  data.frame(from = as.integer(parts[1]),
             to = as.integer(parts[2]),
             weight = as.numeric(parts[3]))
}
edges_df <- do.call(rbind, lapply(edges_raw, parse_edge))

# Create graph
g <- graph_from_data_frame(edges_df, directed = FALSE)

# Identify leaf nodes (0 to n-1)
leaf_nodes <- as.character(0:(n_leaves - 1))

# Compute pairwise distances between leaves
dist_matrix <- matrix(0, nrow = n_leaves, ncol = n_leaves)
for (i in seq_along(leaf_nodes)) {
  dists <- distances(g, v = leaf_nodes[i], to = leaf_nodes, weights = E(g)$weight)
  dist_matrix[i, ] <- dists
}

# Print space-separated matrix
apply(dist_matrix, 1, function(row) cat(paste(row, collapse = "\t"), "\n"))



#### limb length

compute_limb_length <- function(n, j, D) {
  limb_lengths <- c()
  
  for (i in 0:(n - 1)) {
    for (k in 0:(n - 1)) {
      if (i != j && k != j && i != k) {
        limb <- (D[i + 1, j + 1] + D[j + 1, k + 1] - D[i + 1, k + 1]) / 2
        limb_lengths <- c(limb_lengths, limb)
      }
    }
  }
  
  min_limb <- min(limb_lengths)
  return(min_limb)
}

# Sample input
n <- 4
j <- 1
D <- ""
D <- unlist(strsplit(D," |\t|\n"))
D <- as.integer(D)
D <- matrix(D, nrow = n, byrow = TRUE)
D <- matrix(
  c(0, 14, 17, 17,
    14, 0, 7, 13,
    17, 7, 0, 16,
    17, 13, 16, 0),
  nrow = n,
  byrow = TRUE
)
# Compute limb length for leaf j
limb_length <- compute_limb_length_On(j, D)
cat("Limb length for leaf", j, "is:", limb_length, "\n")
# verified

compute_limb_length_On <- function(j, D, k = NULL) {
  n <- nrow(D)
  stopifnot(n == ncol(D))
  stopifnot(j >= 1 && j <= n)
  
  # pick any k != j if not provided
  if (is.null(k)) {
    k <- if (j != 1) 1 else 2
  }
  if (k == j) {
    k <- if (j != 1) 1 else 2
  }
  
  i_idx <- setdiff(seq_len(n), j)
  L <- (D[i_idx, j] + D[j, k] - D[i_idx, k]) / 2
  
  limb <- min(L)
  i_star <- i_idx[which.min(L)]
  
  list(limb = limb, sibling = i_star, k_used = k, L_by_i = setNames(L, i_idx))
}