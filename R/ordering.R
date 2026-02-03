#' Topological Ordering Selection
#'
#' @description
#' Functions for enumerating and selecting topological orderings
#' for the multilinear inequality test.
#'
#' @name ordering
NULL

#' Get All Topological Orderings
#'
#' @description
#' Enumerates all valid topological orderings of a DAG.
#'
#' @param dag A dag object
#' @param max_orderings Maximum number of orderings to return (default 100)
#'
#' @return A list of character vectors, each representing a valid topological
#'   ordering.
#'
#' @details
#' A topological ordering of a DAG is a linear ordering of vertices such that

#' for every directed edge (u, v), vertex u comes before v in the ordering.
#'
#' The number of valid orderings can grow factorially with the number of nodes.
#' Use \code{max_orderings} to limit the enumeration.
#'
#' @examples
#' # Chain: only one ordering
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#' get_topological_orders(g)
#'
#' # Fork: multiple orderings
#' A <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#' get_topological_orders(g)
#'
#' @export
get_topological_orders <- function(dag, max_orderings = 100) {
  # This is an alias for topological_orders in dag_class.R
  topological_orders(dag, max_orderings = max_orderings)
}

#' Select Optimal Ordering
#'
#' @description
#' Selects a topological ordering for the test based on data-driven criteria.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object
#' @param method Selection method:
#'   \itemize{
#'     \item "optimal": Select ordering with maximum Q_n (default)
#'     \item "first": Use first valid ordering
#'     \item "random": Random valid ordering
#'     \item "bonferroni": Test all orderings with Bonferroni correction
#'     \item "fisher": Combine p-values using Fisher's method
#'   }
#' @param max_orderings Maximum orderings to consider for "optimal" and "bonferroni"
#' @param ... Additional arguments passed to compute_qn
#'
#' @return A list with:
#'   \item{ordering}{The selected ordering}
#'   \item{method}{Method used}
#'   \item{qn}{Q_n value for selected ordering (for "optimal")}
#'   \item{all_qn}{Q_n values for all orderings considered (for "optimal")}
#'
#' @details
#' The choice of topological ordering can affect the test's power. The "optimal"
#' method selects the ordering that maximizes Q_n, which tends to give higher
#' power against alternatives.
#'
#' For testing multiple orderings, "bonferroni" applies a Bonferroni correction
#' to control the family-wise error rate, while "fisher" combines individual
#' p-values using Fisher's combination method.
#'
#' @examples
#' A <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' set.seed(123)
#' data <- cbind(rnorm(100), rnorm(100), rnorm(100))
#'
#' select_ordering(data, g, method = "optimal")
#'
#' @export
select_ordering <- function(data, dag, method = "optimal",
                             max_orderings = 100, ...) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  method <- match.arg(method, c("optimal", "first", "random",
                                 "bonferroni", "fisher"))

  # Get all orderings
  orderings <- get_topological_orders(dag, max_orderings = max_orderings)
  n_orderings <- length(orderings)

  if (n_orderings == 0) {
    stop("No valid topological orderings found")
  }

  if (method == "first") {
    return(list(
      ordering = orderings[[1]],
      method = method,
      n_orderings = n_orderings
    ))
  }

  if (method == "random") {
    idx <- sample(n_orderings, 1)
    return(list(
      ordering = orderings[[idx]],
      method = method,
      n_orderings = n_orderings
    ))
  }

  if (method == "optimal") {
    # Compute Q_n for each ordering
    qn_values <- numeric(n_orderings)

    for (i in seq_len(n_orderings)) {
      qn_result <- compute_qn(data, dag, ordering = orderings[[i]], ...)
      qn_values[i] <- qn_result$qn
    }

    best_idx <- which.max(qn_values)
    names(qn_values) <- sapply(orderings, paste, collapse = "->")

    return(list(
      ordering = orderings[[best_idx]],
      method = method,
      qn = qn_values[best_idx],
      all_qn = qn_values,
      n_orderings = n_orderings
    ))
  }

  # For bonferroni and fisher, return all orderings with info
  # Actual testing happens in mi_test
  list(
    orderings = orderings,
    method = method,
    n_orderings = n_orderings
  )
}

#' Compute Number of Topological Orderings
#'
#' @description
#' Computes or estimates the number of valid topological orderings for a DAG.
#'
#' @param dag A dag object
#' @param exact Logical; compute exact count (can be slow) or estimate
#'
#' @return Integer count of topological orderings
#'
#' @details
#' Computing the exact number of topological orderings is #P-complete,
#' so for large DAGs an estimate may be more practical.
#'
#' @export
count_orderings <- function(dag, exact = TRUE) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  if (exact) {
    # Use dynamic programming for small DAGs
    if (dag$n <= 15) {
      return(.count_orderings_exact(dag))
    } else {
      warning("DAG too large for exact count; using estimation")
      exact <- FALSE
    }
  }

  if (!exact) {
    # Estimate using sampling
    return(.estimate_ordering_count(dag))
  }
}

#' Exact Ordering Count via DP
#' @keywords internal
.count_orderings_exact <- function(dag) {
  n <- dag$n
  adj <- dag$adjacency

  # Dynamic programming on subsets
  # dp[S] = number of orderings of subset S
  dp <- rep(0, 2^n)
  dp[1] <- 1  # Empty set

  for (mask in 0:(2^n - 1)) {
    if (dp[mask + 1] == 0) next

    # Find nodes that can be added (no parent in remaining)
    remaining <- which(bitwAnd(mask, 2^(0:(n-1))) == 0)

    for (node in remaining) {
      # Check if all parents are in current set
      parents_in <- which(adj[, node] != 0)
      all_parents_included <- all(bitwAnd(mask, 2^(parents_in - 1)) > 0)

      if (length(parents_in) == 0 || all_parents_included) {
        new_mask <- bitwOr(mask, 2^(node - 1))
        dp[new_mask + 1] <- dp[new_mask + 1] + dp[mask + 1]
      }
    }
  }

  dp[2^n]
}

#' Estimate Ordering Count
#' @keywords internal
.estimate_ordering_count <- function(dag, n_samples = 1000) {
  # Use uniform random ordering generation and counting
  # This is a rough estimate

  n <- dag$n

  # Generate many random orderings
  count <- 0
  for (i in seq_len(n_samples)) {
    # Random permutation
    perm <- sample(dag$nodes)
    if (is_valid_ordering(dag, perm)) {
      count <- count + 1
    }
  }

  # Estimate: (count/n_samples) * n!
  round((count / n_samples) * factorial(n))
}

#' Bonferroni Correction for Multiple Orderings
#'
#' @description
#' Applies Bonferroni correction when testing multiple orderings.
#'
#' @param p_values Vector of p-values from different orderings
#' @param alpha Significance level (default 0.05)
#'
#' @return A list with:
#'   \item{adjusted_p}{Adjusted p-values}
#'   \item{min_p}{Minimum p-value}
#'   \item{reject}{Which hypotheses are rejected}
#'   \item{global_reject}{Whether to reject globally}
#'
#' @export
bonferroni_orderings <- function(p_values, alpha = 0.05) {
  n <- length(p_values)
  adjusted_p <- pmin(p_values * n, 1)
  min_p <- min(p_values)

  list(
    adjusted_p = adjusted_p,
    min_p = min_p,
    adjusted_min_p = min(min_p * n, 1),
    reject = adjusted_p < alpha,
    global_reject = min_p < alpha / n,
    alpha = alpha,
    n_tests = n
  )
}

#' Fisher's Combination Method
#'
#' @description
#' Combines p-values from multiple orderings using Fisher's method.
#'
#' @param p_values Vector of p-values
#'
#' @return A list with:
#'   \item{combined_p}{Combined p-value}
#'   \item{chi_sq}{Chi-squared statistic}
#'   \item{df}{Degrees of freedom}
#'
#' @details
#' Fisher's method combines p-values using:
#' X = -2 * sum(log(p_i))
#' which follows a chi-squared distribution with 2k degrees of freedom
#' under the null hypothesis of all nulls being true.
#'
#' @export
fisher_combine <- function(p_values) {
  # Handle edge cases
  p_values[p_values < 1e-15] <- 1e-15
  p_values[p_values > 1 - 1e-15] <- 1 - 1e-15

  k <- length(p_values)
  chi_sq <- -2 * sum(log(p_values))
  df <- 2 * k

  combined_p <- pchisq(chi_sq, df, lower.tail = FALSE)

  list(
    combined_p = combined_p,
    chi_sq = chi_sq,
    df = df,
    n_tests = k
  )
}

#' Stouffer's Combination Method
#'
#' @description
#' Combines p-values using Stouffer's Z-score method.
#'
#' @param p_values Vector of p-values
#' @param weights Optional weights for each p-value
#'
#' @return A list with:
#'   \item{combined_p}{Combined p-value}
#'   \item{z_score}{Combined Z-score}
#'
#' @export
stouffer_combine <- function(p_values, weights = NULL) {
  # Handle edge cases
  p_values[p_values < 1e-15] <- 1e-15
  p_values[p_values > 1 - 1e-15] <- 1 - 1e-15

  k <- length(p_values)

  if (is.null(weights)) {
    weights <- rep(1, k)
  }

  # Convert to Z-scores
  z_scores <- qnorm(1 - p_values)

  # Weighted combination
  z_combined <- sum(weights * z_scores) / sqrt(sum(weights^2))

  combined_p <- pnorm(z_combined, lower.tail = FALSE)

  list(
    combined_p = combined_p,
    z_score = z_combined,
    n_tests = k
  )
}
