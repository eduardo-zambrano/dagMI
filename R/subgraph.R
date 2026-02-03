#' Subgraph Testing
#'
#' @description
#' Functions for extracting and testing subgraphs of a DAG.
#'
#' @name subgraph
NULL

#' Extract Subgraph
#'
#' @description
#' Extracts a subgraph containing specified nodes from a DAG.
#'
#' @param dag A dag object
#' @param nodes Character vector of node names to include
#'
#' @return A dag object representing the subgraph
#'
#' @details
#' The subgraph contains only the specified nodes and edges between them.
#' Nodes that were connected through intermediate nodes in the original
#' DAG will NOT have direct edges in the subgraph.
#'
#' @examples
#' # Create DAG: X1 -> X2 -> X3 -> X4
#' A <- matrix(c(0, 1, 0, 0,
#'               0, 0, 1, 0,
#'               0, 0, 0, 1,
#'               0, 0, 0, 0), 4, 4, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))
#'
#' # Extract subgraph with X1, X2, X3
#' sub_g <- extract_subgraph(g, c("X1", "X2", "X3"))
#'
#' @export
extract_subgraph <- function(dag, nodes) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  # Validate nodes
  if (!all(nodes %in% dag$nodes)) {
    missing <- setdiff(nodes, dag$nodes)
    stop("Nodes not in DAG: ", paste(missing, collapse = ", "))
  }

  # Get indices
  idx <- match(nodes, dag$nodes)

  # Extract submatrix
  sub_adj <- dag$adjacency[idx, idx, drop = FALSE]
  rownames(sub_adj) <- colnames(sub_adj) <- nodes

  dag(sub_adj, nodes = nodes)
}

#' Extract Ancestral Subgraph
#'
#' @description
#' Extracts a subgraph containing specified nodes and all their ancestors.
#'
#' @param dag A dag object
#' @param nodes Character vector of target node names
#'
#' @return A dag object representing the ancestral subgraph
#'
#' @examples
#' A <- matrix(c(0, 1, 0, 0,
#'               0, 0, 1, 0,
#'               0, 0, 0, 1,
#'               0, 0, 0, 0), 4, 4, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))
#'
#' # Get ancestral subgraph of X4
#' anc_g <- extract_ancestral_subgraph(g, "X4")  # Includes all nodes
#'
#' @export
extract_ancestral_subgraph <- function(dag, nodes) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  # Get all ancestors
  all_nodes <- unique(nodes)
  for (node in nodes) {
    all_nodes <- unique(c(all_nodes, ancestors(dag, node)))
  }

  # Preserve original ordering
  all_nodes <- dag$nodes[dag$nodes %in% all_nodes]

  extract_subgraph(dag, all_nodes)
}

#' Test Multiple Subgraphs
#'
#' @description
#' Tests multiple subgraphs of a DAG with appropriate multiple testing
#' correction.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object
#' @param subgraph_nodes List of character vectors, each specifying nodes
#'   for a subgraph. If NULL, uses all k-node subgraphs.
#' @param k Size of subgraphs to test (used if subgraph_nodes is NULL)
#' @param correction Multiple testing correction: "bonferroni" (default),
#'   "holm", "bh", or "none"
#' @param B Number of bootstrap replicates
#' @param alpha Significance level
#' @param parallel Use parallel processing
#' @param verbose Show progress
#'
#' @return An object of class "subgraph_test_result" with:
#'   \item{results}{List of mi_test_result objects for each subgraph}
#'   \item{p_values}{Vector of raw p-values}
#'   \item{adjusted_p}{Vector of adjusted p-values}
#'   \item{rejected}{Logical vector indicating which subgraphs are rejected}
#'   \item{correction}{Correction method used}
#'
#' @examples
#' A <- matrix(c(0, 1, 0, 0,
#'               0, 0, 1, 0,
#'               0, 0, 0, 1,
#'               0, 0, 0, 0), 4, 4, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))
#'
#' set.seed(123)
#' data <- cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100))
#'
#' # Test all 3-node subgraphs
#' result <- subgraph_test(data, g, k = 3, B = 50)
#'
#' @export
subgraph_test <- function(data, dag, subgraph_nodes = NULL, k = 3,
                           correction = "bonferroni",
                           B = 500, alpha = 0.05,
                           parallel = FALSE, verbose = TRUE) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  correction <- match.arg(correction, c("bonferroni", "holm", "bh", "none"))

  # Generate subgraph node sets if not provided
  if (is.null(subgraph_nodes)) {
    subgraph_nodes <- .generate_subgraphs(dag, k)
  }

  n_subgraphs <- length(subgraph_nodes)

  if (verbose) {
    cat("Testing", n_subgraphs, "subgraphs\n")
  }

  # Test each subgraph
  results <- vector("list", n_subgraphs)
  p_values <- numeric(n_subgraphs)

  for (i in seq_len(n_subgraphs)) {
    nodes_i <- subgraph_nodes[[i]]

    if (verbose) {
      cat("  Subgraph", i, ":", paste(nodes_i, collapse = ", "), "\n")
    }

    # Extract subgraph
    sub_dag <- extract_subgraph(dag, nodes_i)

    # Extract relevant data columns
    sub_data <- data[, nodes_i, drop = FALSE]

    # Test
    results[[i]] <- mi_test(
      sub_data, sub_dag,
      B = B, alpha = alpha,
      parallel = parallel,
      verbose = FALSE
    )

    p_values[i] <- results[[i]]$p_value
  }

  # Apply multiple testing correction
  adjusted_p <- switch(correction,
                        bonferroni = pmin(p_values * n_subgraphs, 1),
                        holm = .holm_correction(p_values),
                        bh = .bh_correction(p_values),
                        none = p_values
  )

  rejected <- adjusted_p < alpha

  if (verbose) {
    cat("\nResults (", correction, " correction):\n", sep = "")
    for (i in seq_len(n_subgraphs)) {
      status <- if (rejected[i]) "REJECTED" else "not rejected"
      cat("  ", paste(subgraph_nodes[[i]], collapse = ","), ": p = ",
          format(p_values[i], digits = 3), " -> ", status, "\n", sep = "")
    }
  }

  structure(
    list(
      results = results,
      subgraph_nodes = subgraph_nodes,
      p_values = p_values,
      adjusted_p = adjusted_p,
      rejected = rejected,
      correction = correction,
      alpha = alpha,
      n_subgraphs = n_subgraphs
    ),
    class = "subgraph_test_result"
  )
}

#' Generate All k-Node Subgraphs
#' @keywords internal
.generate_subgraphs <- function(dag, k) {
  if (k > dag$n) {
    stop("k cannot exceed number of nodes in DAG")
  }

  # All combinations of k nodes
  combs <- combn(dag$nodes, k, simplify = FALSE)

  # Filter to connected subgraphs (optional - could keep all)
  # For now, keep all combinations
  combs
}

#' Holm Correction
#' @keywords internal
.holm_correction <- function(p_values) {
  n <- length(p_values)
  ord <- order(p_values)
  adjusted <- numeric(n)

  for (i in seq_len(n)) {
    adjusted[ord[i]] <- min(p_values[ord[i]] * (n - i + 1), 1)
  }

  # Enforce monotonicity
  for (i in 2:n) {
    adjusted[ord[i]] <- max(adjusted[ord[i]], adjusted[ord[i - 1]])
  }

  adjusted
}

#' Benjamini-Hochberg Correction
#' @keywords internal
.bh_correction <- function(p_values) {
  n <- length(p_values)
  ord <- order(p_values)
  adjusted <- numeric(n)

  for (i in seq_len(n)) {
    adjusted[ord[i]] <- p_values[ord[i]] * n / i
  }

  # Enforce monotonicity (in reverse)
  for (i in (n - 1):1) {
    adjusted[ord[i]] <- min(adjusted[ord[i]], adjusted[ord[i + 1]])
  }

  pmin(adjusted, 1)
}

#' Print Method for Subgraph Test
#' @export
print.subgraph_test_result <- function(x, ...) {
  cat("\n")
  cat("Subgraph Testing Results\n")
  cat("========================\n\n")

  cat("Number of subgraphs tested:", x$n_subgraphs, "\n")
  cat("Correction method:", x$correction, "\n")
  cat("Significance level:", x$alpha, "\n\n")

  cat("Results:\n")
  for (i in seq_len(x$n_subgraphs)) {
    status <- if (x$rejected[i]) "*REJECTED*" else "not rejected"
    cat("  ", paste(x$subgraph_nodes[[i]], collapse = ", "), "\n",
        "    raw p =", format(x$p_values[i], digits = 4),
        ", adj p =", format(x$adjusted_p[i], digits = 4),
        " [", status, "]\n", sep = "")
  }

  cat("\n")
  cat("Rejected:", sum(x$rejected), "/", x$n_subgraphs, "subgraphs\n")

  invisible(x)
}

#' Find Problematic Subgraphs
#'
#' @description
#' Identifies subgraphs that show evidence against DAG compatibility.
#'
#' @param subgraph_result A subgraph_test_result object
#'
#' @return List of subgraph node sets that were rejected
#'
#' @export
problematic_subgraphs <- function(subgraph_result) {
  if (!inherits(subgraph_result, "subgraph_test_result")) {
    stop("Input must be a subgraph_test_result object")
  }

  idx <- which(subgraph_result$rejected)

  if (length(idx) == 0) {
    message("No subgraphs rejected")
    return(list())
  }

  subgraph_result$subgraph_nodes[idx]
}

#' Identify Common Nodes in Rejected Subgraphs
#'
#' @description
#' Finds nodes that appear frequently in rejected subgraphs, which may
#' indicate problem areas in the DAG specification.
#'
#' @param subgraph_result A subgraph_test_result object
#'
#' @return Data frame with nodes and their frequency in rejected subgraphs
#'
#' @export
problematic_nodes <- function(subgraph_result) {
  rejected_subs <- problematic_subgraphs(subgraph_result)

  if (length(rejected_subs) == 0) {
    return(data.frame(node = character(), count = integer(),
                       proportion = numeric()))
  }

  all_nodes <- unlist(rejected_subs)
  node_counts <- table(all_nodes)

  data.frame(
    node = names(node_counts),
    count = as.integer(node_counts),
    proportion = as.numeric(node_counts) / length(rejected_subs),
    stringsAsFactors = FALSE
  )[order(-node_counts), ]
}
