#' Create a DAG Object
#'
#' @description
#' Creates an S3 object representing a directed acyclic graph (DAG) from
#' an adjacency matrix or edge list.
#'
#' @param adjacency A square adjacency matrix where \code{A[i,j] = 1} indicates
#'   an edge from node i to node j, or NULL if using \code{edges}.
#' @param edges A two-column matrix or data frame of edges, where each row
#'   represents an edge from the first column to the second column.
#'   Can be node names or indices.
#' @param nodes Character vector of node names. If NULL, uses column names of
#'   adjacency matrix or generates names X1, X2, ...
#'
#' @return An S3 object of class "dag" with components:
#'   \item{adjacency}{The adjacency matrix}
#'   \item{nodes}{Character vector of node names}
#'   \item{n}{Number of nodes}
#'
#' @details
#' The adjacency matrix convention is that \code{A[i,j] = 1} means there is
#' a directed edge from node i to node j (i.e., i is a parent of j).
#' The function validates that the graph is acyclic.
#'
#' @examples
#' # Create a chain DAG: X1 -> X2 -> X3
#' A <- matrix(c(0, 1, 0,
#'               0, 0, 1,
#'               0, 0, 0), nrow = 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' # Create from edge list
#' edges <- data.frame(from = c("X1", "X2"), to = c("X2", "X3"))
#' g <- dag(edges = edges)
#'
#' @export
dag <- function(adjacency = NULL, edges = NULL, nodes = NULL) {
  # Validate inputs

if (is.null(adjacency) && is.null(edges)) {
    stop("Must provide either 'adjacency' matrix or 'edges'")
  }

  if (!is.null(adjacency) && !is.null(edges)) {
    stop("Provide only one of 'adjacency' or 'edges', not both")
  }

  # Build from edge list if provided
  if (!is.null(edges)) {
    edges <- as.data.frame(edges)
    if (ncol(edges) < 2) {
      stop("'edges' must have at least 2 columns (from, to)")
    }

    # Get unique nodes
    all_nodes <- unique(c(as.character(edges[[1]]), as.character(edges[[2]])))
    if (is.null(nodes)) {
      nodes <- sort(all_nodes)
    } else {
      # Ensure all edge nodes are in nodes
      missing <- setdiff(all_nodes, nodes)
      if (length(missing) > 0) {
        stop("Edges reference nodes not in 'nodes': ", paste(missing, collapse = ", "))
      }
    }

    n <- length(nodes)
    adjacency <- matrix(0, nrow = n, ncol = n)
    rownames(adjacency) <- colnames(adjacency) <- nodes

    for (i in seq_len(nrow(edges))) {
      from_idx <- match(as.character(edges[i, 1]), nodes)
      to_idx <- match(as.character(edges[i, 2]), nodes)
      adjacency[from_idx, to_idx] <- 1
    }
  }

  # Validate adjacency matrix
  adjacency <- as.matrix(adjacency)
  if (!is.numeric(adjacency)) {
    stop("Adjacency matrix must be numeric")
  }
  if (nrow(adjacency) != ncol(adjacency)) {
    stop("Adjacency matrix must be square")
  }

  n <- nrow(adjacency)

  # Set node names
  if (is.null(nodes)) {
    if (!is.null(colnames(adjacency))) {
      nodes <- colnames(adjacency)
    } else {
      nodes <- paste0("X", seq_len(n))
    }
  }

  if (length(nodes) != n) {
    stop("Length of 'nodes' must match dimensions of adjacency matrix")
  }

  rownames(adjacency) <- colnames(adjacency) <- nodes

  # Ensure binary
  adjacency <- (adjacency != 0) * 1

  # Check for self-loops
  if (any(diag(adjacency) != 0)) {
    stop("DAG cannot have self-loops")
  }

  # Check acyclicity using topological sort
  if (!.is_acyclic(adjacency)) {
    stop("Graph contains cycles and is not a DAG")
  }

  structure(
    list(
      adjacency = adjacency,
      nodes = nodes,
      n = n
    ),
    class = "dag"
  )
}

#' Check if Graph is Acyclic
#' @keywords internal
.is_acyclic <- function(adj) {
  n <- nrow(adj)
  if (n == 0) return(TRUE)

  # Kahn's algorithm
  in_degree <- colSums(adj)
  queue <- which(in_degree == 0)
  count <- 0

  adj_copy <- adj

  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]
    count <- count + 1

    # Remove outgoing edges
    children <- which(adj_copy[node, ] != 0)
    for (child in children) {
      adj_copy[node, child] <- 0
      in_degree[child] <- in_degree[child] - 1
      if (in_degree[child] == 0) {
        queue <- c(queue, child)
      }
    }
  }

  count == n
}

#' Test if Object is a DAG
#'
#' @param x An R object
#' @return Logical indicating if x is a dag object
#' @export
is_dag <- function(x) {
  inherits(x, "dag")
}

#' Print Method for DAG
#' @param x A dag object
#' @param ... Additional arguments (ignored)
#' @export
print.dag <- function(x, ...) {
  cat("DAG with", x$n, "nodes:", paste(x$nodes, collapse = ", "), "\n")
  cat("Edges:\n")

  edges <- which(x$adjacency != 0, arr.ind = TRUE)
  if (nrow(edges) == 0) {
    cat("  (no edges)\n")
  } else {
    for (i in seq_len(nrow(edges))) {
      cat("  ", x$nodes[edges[i, 1]], "->", x$nodes[edges[i, 2]], "\n")
    }
  }
  invisible(x)
}

#' Plot Method for DAG
#'
#' @param x A dag object
#' @param layout Layout algorithm: "hierarchical" (default) or "circular"
#' @param ... Additional arguments passed to plotting functions
#' @export
plot.dag <- function(x, layout = "hierarchical", ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    # Simple base R plot
    .plot_dag_base(x, layout, ...)
  } else {
    .plot_dag_igraph(x, layout, ...)
  }
}

#' Base R Plot for DAG
#' @keywords internal
.plot_dag_base <- function(dag, layout = "hierarchical", ...) {
  n <- dag$n

  if (layout == "circular") {
    angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
    x <- cos(angles)
    y <- sin(angles)
  } else {
    # Hierarchical layout based on topological order
    order <- .topological_sort(dag$adjacency)
    levels <- .compute_levels(dag$adjacency, order)

    # Compute x, y positions
    max_level <- max(levels)
    x <- numeric(n)
    y <- (max_level - levels) / max(max_level, 1)

    # Spread nodes at same level
    for (lev in unique(levels)) {
      nodes_at_level <- which(levels == lev)
      if (length(nodes_at_level) > 1) {
        x[nodes_at_level] <- seq(-1, 1, length.out = length(nodes_at_level))
      }
    }
  }

  # Plot
  oldpar <- par(mar = c(1, 1, 1, 1))
  on.exit(par(oldpar))

  plot(x, y, type = "n", xlim = c(-1.5, 1.5), ylim = c(-0.5, 1.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", asp = 1, ...)

  # Draw edges with arrows
  edges <- which(dag$adjacency != 0, arr.ind = TRUE)
  for (i in seq_len(nrow(edges))) {
    from <- edges[i, 1]
    to <- edges[i, 2]
    arrows(x[from], y[from], x[to], y[to],
           length = 0.1, angle = 20, col = "gray40")
  }

  # Draw nodes
  points(x, y, pch = 21, bg = "lightblue", cex = 3)
  text(x, y, dag$nodes, cex = 0.8)
}

#' igraph Plot for DAG
#' @keywords internal
.plot_dag_igraph <- function(dag, layout = "hierarchical", ...) {
  g <- igraph::graph_from_adjacency_matrix(dag$adjacency, mode = "directed")

  if (layout == "hierarchical") {
    lay <- igraph::layout_with_sugiyama(g)$layout
  } else {
    lay <- igraph::layout_in_circle(g)
  }

  igraph::plot.igraph(g, layout = lay,
                      vertex.label = dag$nodes,
                      vertex.color = "lightblue",
                      vertex.size = 30,
                      edge.arrow.size = 0.5,
                      ...)
}

#' Get Parents of a Node
#'
#' @param dag A dag object
#' @param node Node name or index
#' @return Character vector of parent node names
#' @export
parents <- function(dag, node) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  idx <- .node_index(dag, node)
  parent_idx <- which(dag$adjacency[, idx] != 0)
  dag$nodes[parent_idx]
}

#' Get Children of a Node
#'
#' @param dag A dag object
#' @param node Node name or index
#' @return Character vector of child node names
#' @export
children <- function(dag, node) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  idx <- .node_index(dag, node)
  child_idx <- which(dag$adjacency[idx, ] != 0)
  dag$nodes[child_idx]
}

#' Get Ancestors of a Node
#'
#' @param dag A dag object
#' @param node Node name or index
#' @return Character vector of ancestor node names
#' @export
ancestors <- function(dag, node) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  idx <- .node_index(dag, node)

  # BFS/DFS backwards
  visited <- logical(dag$n)
  stack <- which(dag$adjacency[, idx] != 0)

  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]

    if (!visited[current]) {
      visited[current] <- TRUE
      # Add parents of current
      new_parents <- which(dag$adjacency[, current] != 0)
      stack <- c(stack, new_parents[!visited[new_parents]])
    }
  }

  dag$nodes[visited]
}

#' Get Descendants of a Node
#'
#' @param dag A dag object
#' @param node Node name or index
#' @return Character vector of descendant node names
#' @export
descendants <- function(dag, node) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  idx <- .node_index(dag, node)

  # BFS/DFS forwards
  visited <- logical(dag$n)
  stack <- which(dag$adjacency[idx, ] != 0)

  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]

    if (!visited[current]) {
      visited[current] <- TRUE
      # Add children of current
      new_children <- which(dag$adjacency[current, ] != 0)
      stack <- c(stack, new_children[!visited[new_children]])
    }
  }

  dag$nodes[visited]
}

#' Get All Topological Orderings
#'
#' @description
#' Enumerates all valid topological orderings of the DAG.
#'
#' @param dag A dag object
#' @param max_orderings Maximum number of orderings to return (default 100).
#'   Set to Inf for all orderings (may be slow for large DAGs).
#'
#' @return A list of character vectors, each representing a valid topological ordering
#'
#' @details
#' For a DAG with n nodes, there can be up to n! orderings. This function
#' uses a recursive algorithm with pruning to enumerate them efficiently.
#'
#' @examples
#' A <- matrix(c(0, 1, 1,
#'               0, 0, 1,
#'               0, 0, 0), nrow = 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#' topological_orders(g)
#'
#' @export
topological_orders <- function(dag, max_orderings = 100) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  orderings <- list()
  .enumerate_orderings(dag$adjacency, dag$nodes,
                       integer(0), orderings, max_orderings)
}

#' Recursive Enumeration of Topological Orderings
#' @keywords internal
.enumerate_orderings <- function(adj, nodes, current, orderings_env, max_ord) {
  n <- nrow(adj)

  if (length(current) == n) {
    # Found a complete ordering
    orderings_env$list <- c(orderings_env$list, list(nodes[current]))
    return(length(orderings_env$list) < max_ord)
  }

  # Find nodes with no incoming edges from remaining nodes
  remaining <- setdiff(seq_len(n), current)
  in_degree <- colSums(adj[remaining, remaining, drop = FALSE])

  # Map back to original indices
  available <- remaining[in_degree == 0]

  for (node in available) {
    if (!.enumerate_orderings(adj, nodes, c(current, node),
                              orderings_env, max_ord)) {
      return(FALSE)
    }
  }

  TRUE
}

# Wrapper to manage environment
topological_orders <- function(dag, max_orderings = 100) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  env <- new.env()
  env$list <- list()

  .enumerate_orderings_env(dag$adjacency, dag$nodes,
                           integer(0), env, max_orderings)

  env$list
}

.enumerate_orderings_env <- function(adj, nodes, current, env, max_ord) {
  n <- nrow(adj)

  if (length(current) == n) {
    env$list <- c(env$list, list(nodes[current]))
    return(length(env$list) < max_ord)
  }

  remaining <- setdiff(seq_len(n), current)

  # Check in-degree from remaining nodes only
  available <- c()
  for (r in remaining) {
    parents_in_remaining <- sum(adj[remaining, r])
    if (parents_in_remaining == 0) {
      available <- c(available, r)
    }
  }

  for (node in available) {
    if (!.enumerate_orderings_env(adj, nodes, c(current, node), env, max_ord)) {
      return(FALSE)
    }
  }

  TRUE
}

#' Check if Ordering is Valid
#'
#' @param dag A dag object
#' @param ordering Character vector of node names representing an ordering
#' @return Logical indicating if the ordering is a valid topological order
#' @export
is_valid_ordering <- function(dag, ordering) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  if (length(ordering) != dag$n) return(FALSE)
  if (!all(ordering %in% dag$nodes)) return(FALSE)
  if (length(unique(ordering)) != length(ordering)) return(FALSE)

  # Check that parents come before children
  for (i in seq_len(dag$n)) {
    node <- ordering[i]
    node_parents <- parents(dag, node)
    if (length(node_parents) > 0) {
      parent_positions <- match(node_parents, ordering)
      if (any(parent_positions >= i)) {
        return(FALSE)
      }
    }
  }

  TRUE
}

#' Convert Node to Index
#' @keywords internal
.node_index <- function(dag, node) {
  if (is.numeric(node)) {
    if (node < 1 || node > dag$n) {
      stop("Node index out of range")
    }
    return(as.integer(node))
  }

  idx <- match(node, dag$nodes)
  if (is.na(idx)) {
    stop("Node '", node, "' not found in DAG")
  }
  idx
}

#' Topological Sort (Single Ordering)
#' @keywords internal
.topological_sort <- function(adj) {
  n <- nrow(adj)
  in_degree <- colSums(adj)
  queue <- which(in_degree == 0)
  order <- integer(n)
  pos <- 1

  adj_copy <- adj

  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]
    order[pos] <- node
    pos <- pos + 1

    children <- which(adj_copy[node, ] != 0)
    for (child in children) {
      adj_copy[node, child] <- 0
      in_degree[child] <- in_degree[child] - 1
      if (in_degree[child] == 0) {
        queue <- c(queue, child)
      }
    }
  }

  order
}

#' Compute Levels for Hierarchical Layout
#' @keywords internal
.compute_levels <- function(adj, order) {
  n <- nrow(adj)
  levels <- integer(n)

  for (node in order) {
    parent_idx <- which(adj[, node] != 0)
    if (length(parent_idx) == 0) {
      levels[node] <- 0
    } else {
      levels[node] <- max(levels[parent_idx]) + 1
    }
  }

  levels
}
