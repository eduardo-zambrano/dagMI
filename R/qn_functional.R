#' Compute Q_n^G Functional
#'
#' @description
#' Computes the DAG-specific functional Q_n^G that appears in Carbery's
#' inequality. This is the key quantity that determines the bound in
#' the multilinear inequality.
#'
#' @param data Numeric matrix with n observations and p variables (columns
#'   ordered according to a topological ordering of the DAG)
#' @param dag A dag object specifying the DAG structure
#' @param ordering Optional topological ordering to use. If NULL, uses the
#'   first valid ordering.
#' @param method Integration method: "quadrature" (default for n <= 5) or
#'   "montecarlo" (default for n > 5)
#' @param n_quad Number of quadrature points per dimension (default 10)
#' @param n_mc Number of Monte Carlo samples (default 10000)
#' @param bandwidth Vector of bandwidths for KDE. If NULL, uses Silverman's rule.
#'
#' @return An object of class "qn_result" with components:
#'   \item{qn}{The computed Q_n^G value}
#'   \item{se}{Standard error (for Monte Carlo method)}
#'   \item{method}{Integration method used}
#'   \item{n_vars}{Number of variables}
#'   \item{ordering}{Topological ordering used}
#'   \item{kde_list}{List of KDE objects (for reuse)}
#'
#' @details
#' Q_n^G is defined as:
#' \deqn{Q_n^G = \int \prod_{j=2}^n f_{j-1,j}(x_{j-1}, x_j)^{1/(n+1)}
#'       \cdot \prod_{j=1}^n f_j(x_j)^{-(n-1)/(n+1)} dx}
#'
#' where f_j is the marginal density of X_j and f_{j-1,j} is the joint
#' density of consecutive pairs in the ordering.
#'
#' For small n (n <= 5), tensor-product Gauss-Legendre quadrature is used.
#' For larger n, importance sampling Monte Carlo is more efficient.
#'
#' @examples
#' # Create simple chain DAG: X1 -> X2 -> X3
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' # Generate data
#' set.seed(123)
#' X1 <- rnorm(200)
#' X2 <- 0.5 * X1 + rnorm(200, sd = 0.5)
#' X3 <- 0.5 * X2 + rnorm(200, sd = 0.5)
#' data <- cbind(X1, X2, X3)
#'
#' # Compute Q_n
#' result <- compute_qn(data, g)
#'
#' @export
compute_qn <- function(data, dag, ordering = NULL,
                        method = NULL, n_quad = 10, n_mc = 10000,
                        bandwidth = NULL) {
  # Validate inputs
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  n_obs <- nrow(data)
  n_vars <- ncol(data)

  # Get ordering
  if (is.null(ordering)) {
    ordering <- topological_orders(dag, max_orderings = 1)[[1]]
  } else {
    if (!is_valid_ordering(dag, ordering)) {
      stop("Invalid topological ordering")
    }
  }

  # Reorder data according to ordering
  data_ordered <- data[, ordering, drop = FALSE]

  # Select method
  if (is.null(method)) {
    method <- if (n_vars <= 5) "quadrature" else "montecarlo"
  }
  method <- match.arg(method, c("quadrature", "montecarlo"))

  # Create KDE objects
  kde_info <- create_kde_list(data_ordered, bandwidth)
  kde_list <- kde_info$kde_objects

  # Compute bounds for integration
  bounds <- matrix(NA, nrow = n_vars, ncol = 2)
  for (j in seq_len(n_vars)) {
    bw_j <- kde_list[[j]]$bandwidth
    bounds[j, ] <- range(data_ordered[, j]) + c(-4, 4) * bw_j
  }

  if (method == "quadrature") {
    result <- .compute_qn_quadrature(data_ordered, kde_list, bounds,
                                      n_vars, n_quad)
  } else {
    result <- .compute_qn_montecarlo(data_ordered, kde_list, bounds,
                                      n_vars, n_mc)
  }

  structure(
    list(
      qn = result$qn,
      se = result$se,
      method = method,
      n_vars = n_vars,
      ordering = ordering,
      kde_list = kde_list,
      bounds = bounds
    ),
    class = "qn_result"
  )
}

#' Quadrature-based Q_n Computation
#' @keywords internal
.compute_qn_quadrature <- function(data, kde_list, bounds, n_vars, n_quad) {
  # Get Gauss-Legendre rule
  gl <- gauss_legendre(n_quad)

  # Transform to each dimension
  rules <- lapply(seq_len(n_vars), function(j) {
    transform_quadrature(gl, bounds[j, 1], bounds[j, 2])
  })

  # Create tensor product grid
  grids <- lapply(rules, function(r) r$nodes)
  wts <- lapply(rules, function(r) r$weights)

  grid <- as.matrix(expand.grid(grids))
  weight_grid <- expand.grid(wts)
  weights <- apply(weight_grid, 1, prod)

  # Call C++ for fast computation
  qn <- compute_qn_quadrature_cpp(kde_list, grid, weights, n_vars)

  list(qn = qn, se = NA)
}

#' Monte Carlo Q_n Computation
#' @keywords internal
.compute_qn_montecarlo <- function(data, kde_list, bounds, n_vars, n_mc) {
  # Generate samples from proposal (uniform over range)
  samples <- matrix(NA, nrow = n_mc, ncol = n_vars)

  for (j in seq_len(n_vars)) {
    samples[, j] <- runif(n_mc, bounds[j, 1], bounds[j, 2])
  }

  # Importance weight: uniform proposal, so weight is volume
  volume <- prod(bounds[, 2] - bounds[, 1])

  # Call C++ for fast computation
  result <- compute_qn_montecarlo_cpp(kde_list, samples, n_vars, n_mc)

  # Multiply by volume (importance sampling correction)
  qn <- result$qn * volume
  se <- result$se * volume

  list(qn = qn, se = se)
}

#' Print Method for qn_result
#' @export
print.qn_result <- function(x, ...) {
  cat("Q_n^G Computation Result\n")
  cat("------------------------\n")
  cat("Q_n value:", format(x$qn, digits = 4), "\n")
  if (!is.na(x$se)) {
    cat("Std. Error:", format(x$se, digits = 4), "\n")
  }
  cat("Method:", x$method, "\n")
  cat("Variables:", x$n_vars, "\n")
  cat("Ordering:", paste(x$ordering, collapse = " -> "), "\n")
  invisible(x)
}

#' Compute Q_n for Multiple Orderings
#'
#' @description
#' Computes Q_n^G for multiple (or all) topological orderings of a DAG.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object
#' @param orderings List of orderings to evaluate, or "all" for all orderings
#' @param max_orderings Maximum number of orderings to evaluate (default 100)
#' @param ... Additional arguments passed to compute_qn
#'
#' @return A list with components:
#'   \item{results}{List of qn_result objects}
#'   \item{qn_values}{Named numeric vector of Q_n values}
#'   \item{best_ordering}{Ordering with maximum Q_n}
#'   \item{best_qn}{Maximum Q_n value}
#'
#' @export
compute_qn_all_orderings <- function(data, dag, orderings = "all",
                                       max_orderings = 100, ...) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  if (identical(orderings, "all")) {
    orderings <- topological_orders(dag, max_orderings = max_orderings)
  }

  n_orderings <- length(orderings)

  results <- vector("list", n_orderings)
  qn_values <- numeric(n_orderings)

  for (i in seq_len(n_orderings)) {
    results[[i]] <- compute_qn(data, dag, ordering = orderings[[i]], ...)
    qn_values[i] <- results[[i]]$qn
  }

  # Name by ordering
  names(qn_values) <- sapply(orderings, paste, collapse = "->")
  names(results) <- names(qn_values)

  # Find best
  best_idx <- which.max(qn_values)

  list(
    results = results,
    qn_values = qn_values,
    best_ordering = orderings[[best_idx]],
    best_qn = qn_values[best_idx]
  )
}

#' Theoretical Bounds on Q_n
#'
#' @description
#' Computes theoretical bounds on Q_n^G based on the DAG structure.
#'
#' @param dag A dag object
#'
#' @return A list with:
#'   \item{lower}{Lower bound (always 0)}
#'   \item{upper}{Upper bound based on DAG structure}
#'   \item{description}{Description of the bounds}
#'
#' @details
#' For a DAG with n nodes:
#' - Q_n >= 0 (always)
#' - Q_n <= 1 when all variables are independent
#' - Q_n = 1 exactly when the data satisfies the DAG constraints with equality
#'
#' @export
qn_bounds <- function(dag) {
  if (!is_dag(dag)) stop("'dag' must be a dag object")

  n <- dag$n

  list(
    lower = 0,
    upper = 1,
    description = paste0(
      "Q_", n, " is bounded in [0, 1] for the multilinear inequality. ",
      "Q_n = 1 occurs when variables satisfy the DAG constraints ",
      "with equality (perfect dependence along edges)."
    )
  )
}
