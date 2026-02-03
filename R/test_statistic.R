#' Test Statistic Computation
#'
#' @description
#' Computes the test statistic T_h^G based on Carbery's inequality with
#' density correction.
#'
#' @name test_statistic
NULL

#' Squared Test Function
#'
#' @description
#' The default test function h(x) = x^2.
#'
#' @param x Numeric vector
#' @return x^2
#'
#' @examples
#' h_squared(c(-1, 0, 1, 2))
#'
#' @export
h_squared <- function(x) {
  x^2
}

#' Polynomial Test Function
#'
#' @description
#' Test function h(x) = x^(2k) for positive integer k.
#'
#' @param x Numeric vector
#' @param k Positive integer power (will be doubled, so h(x) = x^(2k))
#'
#' @return x^(2k)
#'
#' @examples
#' h_poly(c(-1, 0, 1, 2), k = 2)  # x^4
#'
#' @export
h_poly <- function(x, k = 1) {
  if (k < 1 || k != round(k)) {
    stop("'k' must be a positive integer")
  }
  x^(2 * k)
}

#' Exponential Test Function
#'
#' @description
#' Test function h(x) = exp(t * x).
#'
#' @param x Numeric vector
#' @param t Parameter (default 1)
#'
#' @return exp(t * x)
#'
#' @examples
#' h_exp(c(-1, 0, 1, 2), t = 0.5)
#'
#' @export
h_exp <- function(x, t = 1) {
  exp(t * x)
}

#' Indicator Test Function
#'
#' @description
#' Test function h(x) = I(x > c).
#'
#' @param x Numeric vector
#' @param c Threshold (default 0)
#'
#' @return Numeric vector of 0s and 1s
#'
#' @examples
#' h_indicator(c(-1, 0, 1, 2), c = 0.5)
#'
#' @export
h_indicator <- function(x, c = 0) {
  as.numeric(x > c)
}

#' Compute Test Statistic T_h^G
#'
#' @description
#' Computes the test statistic based on Carbery's inequality:
#' T_h^G = RHS - LHS
#' where RHS is the bound from the inequality and LHS is the observed value.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object specifying the DAG structure
#' @param h Test function or name of built-in function ("squared", "poly",
#'   "exp", "indicator"). Default is "squared".
#' @param h_params List of parameters for test function (e.g., list(k = 2) for poly)
#' @param qn_result Optional pre-computed Q_n result. If NULL, computes Q_n.
#' @param ordering Optional topological ordering to use
#' @param bandwidth Vector of bandwidths for KDE. If NULL, uses Silverman's rule.
#'
#' @return An object of class "test_stat_result" with components:
#'   \item{statistic}{The test statistic T_h^G}
#'   \item{lhs}{Left-hand side of inequality (observed)}
#'   \item{rhs}{Right-hand side of inequality (bound)}
#'   \item{qn}{Q_n^G value used}
#'   \item{h_name}{Name of test function}
#'   \item{ordering}{Topological ordering used}
#'   \item{n_obs}{Number of observations}
#'   \item{n_vars}{Number of variables}
#'
#' @details
#' The test statistic is:
#' \deqn{T_h^G = Q_n \cdot \prod_{j=1}^n (E[h_j(X_j)^{n+1}])^{1/(n+1)} -
#'              E\left[\prod_{j=1}^n h_j(X_j) \cdot p_j(X_j)^{1/(n+1)}\right]}
#'
#' Under H_0 (data consistent with DAG), T_h^G >= 0 with high probability.
#' Large negative values suggest violation of the DAG structure.
#'
#' @examples
#' # Create chain DAG and generate data
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' set.seed(123)
#' X1 <- rnorm(200)
#' X2 <- 0.5 * X1 + rnorm(200, sd = 0.5)
#' X3 <- 0.5 * X2 + rnorm(200, sd = 0.5)
#' data <- cbind(X1, X2, X3)
#'
#' # Compute test statistic
#' result <- compute_test_stat(data, g)
#'
#' @export
compute_test_stat <- function(data, dag, h = "squared", h_params = list(),
                               qn_result = NULL, ordering = NULL,
                               bandwidth = NULL) {
  # Validate inputs
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  n_obs <- nrow(data)
  n_vars <- ncol(data)

  # Get/compute Q_n
  if (is.null(qn_result)) {
    qn_result <- compute_qn(data, dag, ordering = ordering,
                             bandwidth = bandwidth)
  }

  qn <- qn_result$qn
  ordering <- qn_result$ordering
  kde_list <- qn_result$kde_list

  # Reorder data
  data_ordered <- data[, ordering, drop = FALSE]

  # Parse test function
  h_func <- .parse_test_function(h, h_params)
  h_name <- if (is.character(h)) h else "custom"

  # Apply test function to each variable
  h_values <- matrix(NA, nrow = n_obs, ncol = n_vars)
  for (j in seq_len(n_vars)) {
    h_values[, j] <- h_func(data_ordered[, j])
  }

  # Compute density correction factor for each observation
  # prod_{j=1}^n p_j(X_j)^{1/(n+1)}
  power_dens <- 1 / (n_vars + 1)
  correction <- compute_density_correction_cpp(data_ordered, kde_list, n_vars)

  # Compute LHS: E[prod h_j(X_j) * prod p_j(X_j)^{1/(n+1)}]
  h_prod <- apply(h_values, 1, prod)
  lhs <- mean(h_prod * correction)

  # Compute RHS: Q_n * prod (E[h_j(X_j)^{n+1}])^{1/(n+1)}
  moments <- compute_moments_cpp(h_values, n_vars, n_obs)
  rhs <- qn * prod(moments)

  # Test statistic
  statistic <- rhs - lhs

  structure(
    list(
      statistic = statistic,
      lhs = lhs,
      rhs = rhs,
      qn = qn,
      h_name = h_name,
      h_values = h_values,
      correction = correction,
      moments = moments,
      ordering = ordering,
      n_obs = n_obs,
      n_vars = n_vars
    ),
    class = "test_stat_result"
  )
}

#' Parse Test Function
#' @keywords internal
.parse_test_function <- function(h, h_params) {
  if (is.function(h)) {
    return(h)
  }

  if (!is.character(h)) {
    stop("'h' must be a function or character string")
  }

  h <- match.arg(h, c("squared", "poly", "exp", "indicator"))

  switch(h,
         "squared" = h_squared,
         "poly" = function(x) do.call(h_poly, c(list(x = x), h_params)),
         "exp" = function(x) do.call(h_exp, c(list(x = x), h_params)),
         "indicator" = function(x) do.call(h_indicator, c(list(x = x), h_params))
  )
}

#' Print Method for test_stat_result
#' @export
print.test_stat_result <- function(x, ...) {
  cat("Test Statistic T_h^G\n")
  cat("--------------------\n")
  cat("T_h^G =", format(x$statistic, digits = 4), "\n")
  cat("  LHS (observed):", format(x$lhs, digits = 4), "\n")
  cat("  RHS (bound):", format(x$rhs, digits = 4), "\n")
  cat("  Q_n:", format(x$qn, digits = 4), "\n")
  cat("Test function:", x$h_name, "\n")
  cat("Ordering:", paste(x$ordering, collapse = " -> "), "\n")
  cat("n =", x$n_obs, ", p =", x$n_vars, "\n")
  invisible(x)
}

#' Compute Test Statistics for Multiple Test Functions
#'
#' @description
#' Computes T_h^G for multiple test functions to assess robustness.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object
#' @param test_functions Named list of test function specifications
#' @param qn_result Optional pre-computed Q_n result
#' @param ordering Optional topological ordering
#'
#' @return A list with:
#'   \item{results}{List of test_stat_result objects}
#'   \item{statistics}{Named vector of test statistics}
#'   \item{summary}{Summary data frame}
#'
#' @examples
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' set.seed(123)
#' data <- cbind(rnorm(200), rnorm(200), rnorm(200))
#'
#' # Multiple test functions
#' funcs <- list(
#'   squared = list(h = "squared"),
#'   poly4 = list(h = "poly", h_params = list(k = 2)),
#'   exp05 = list(h = "exp", h_params = list(t = 0.5))
#' )
#' results <- compute_test_stat_multiple(data, g, funcs)
#'
#' @export
compute_test_stat_multiple <- function(data, dag, test_functions = NULL,
                                         qn_result = NULL, ordering = NULL) {
  # Default test functions
  if (is.null(test_functions)) {
    test_functions <- list(
      squared = list(h = "squared"),
      poly4 = list(h = "poly", h_params = list(k = 2)),
      exp05 = list(h = "exp", h_params = list(t = 0.5))
    )
  }

  # Compute Q_n once
  if (is.null(qn_result)) {
    qn_result <- compute_qn(data, dag, ordering = ordering)
  }

  n_funcs <- length(test_functions)
  results <- vector("list", n_funcs)
  statistics <- numeric(n_funcs)

  for (i in seq_len(n_funcs)) {
    func_spec <- test_functions[[i]]
    results[[i]] <- compute_test_stat(
      data, dag,
      h = func_spec$h,
      h_params = func_spec$h_params %||% list(),
      qn_result = qn_result
    )
    statistics[i] <- results[[i]]$statistic
  }

  names(results) <- names(test_functions)
  names(statistics) <- names(test_functions)

  # Summary data frame
  summary_df <- data.frame(
    test_function = names(test_functions),
    statistic = statistics,
    lhs = sapply(results, function(r) r$lhs),
    rhs = sapply(results, function(r) r$rhs),
    stringsAsFactors = FALSE
  )

  list(
    results = results,
    statistics = statistics,
    summary = summary_df,
    qn = qn_result$qn
  )
}

#' Check Admissibility of Test Function
#'
#' @description
#' Verifies that a test function satisfies the admissibility conditions
#' for the multilinear inequality test.
#'
#' @param h Test function
#' @param data Data matrix to evaluate on
#'
#' @return A list with:
#'   \item{admissible}{Logical indicating if function is admissible}
#'   \item{non_negative}{Logical indicating if h(x) >= 0}
#'   \item{finite_moments}{Logical indicating if moments are finite}
#'   \item{message}{Description of any issues}
#'
#' @details
#' A test function h is admissible if:
#' \enumerate{
#'   \item \code{h(x) >= 0} for all x (non-negativity)
#'   \item \code{E[h(X)^(n+1)] < Inf} (finite moments)
#' }
#'
#' @export
check_admissibility <- function(h, data) {
  data <- as.matrix(data)
  n_vars <- ncol(data)

  # Apply to each column
  h_values <- apply(data, 2, h)

  # Check non-negativity
  non_negative <- all(h_values >= 0)

  # Check finite moments
  moments_finite <- all(is.finite(h_values^(n_vars + 1)))
  mean_finite <- all(is.finite(colMeans(h_values^(n_vars + 1))))

  admissible <- non_negative && moments_finite && mean_finite

  message <- if (admissible) {
    "Test function is admissible"
  } else {
    issues <- c()
    if (!non_negative) issues <- c(issues, "h(x) < 0 for some x")
    if (!moments_finite) issues <- c(issues, "h(x)^(n+1) not finite")
    if (!mean_finite) issues <- c(issues, "E[h(X)^(n+1)] not finite")
    paste("Issues:", paste(issues, collapse = "; "))
  }

  list(
    admissible = admissible,
    non_negative = non_negative,
    finite_moments = moments_finite && mean_finite,
    message = message
  )
}
