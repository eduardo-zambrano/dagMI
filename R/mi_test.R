#' Multilinear Inequality Test for DAG Structures
#'
#' @description
#' Main interface for the multilinear inequality test based on Carbery's
#' inequality.
#'
#' @name mi_test
NULL

#' Test DAG Compatibility Using Multilinear Inequality
#'
#' @description
#' Tests whether observed data is compatible with a specified DAG structure
#' using Carbery's multilinear inequality.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object specifying the DAG structure to test
#' @param test_functions Test function(s) to use. Can be:
#'   \itemize{
#'     \item A character string: "squared" (default), "poly", "exp", "indicator"
#'     \item A function
#'     \item A list of functions or specifications
#'   }
#' @param ordering Topological ordering to use. Can be:
#'   \itemize{
#'     \item "optimal": Select ordering with maximum Q_n (default)
#'     \item "first": Use first valid ordering
#'     \item "all": Test all orderings with multiple testing correction
#'     \item A character vector specifying a particular ordering
#'   }
#' @param B Number of bootstrap replicates (default 500)
#' @param alpha Significance level (default 0.05)
#' @param parallel Logical or number of cores for parallel processing
#' @param verbose Logical; show progress (default TRUE)
#' @param seed Random seed for reproducibility
#'
#' @return An S3 object of class "mi_test_result" with components:
#'   \item{statistic}{Observed test statistic T_h^G}
#'   \item{p_value}{Bootstrap p-value}
#'   \item{decision}{Test decision: "reject" or "fail_to_reject"}
#'   \item{alpha}{Significance level used}
#'   \item{qn}{Q_n^G value}
#'   \item{bootstrap}{Bootstrap results (if B > 0)}
#'   \item{ordering}{Topological ordering used}
#'   \item{dag}{The DAG tested}
#'   \item{n_obs}{Number of observations}
#'   \item{n_vars}{Number of variables}
#'
#' @details
#' The test is based on Carbery's multilinear inequality, which states that
#' for a valid DAG:
#' \deqn{E\left[\prod_{j=1}^n h_j(X_j) \cdot p_j(X_j)^{1/(n+1)}\right] \leq
#'       Q_n^G \prod_{j=1}^n \left(E[h_j(X_j)^{n+1}]\right)^{1/(n+1)}}
#'
#' The test statistic T_h^G = RHS - LHS should be non-negative under H_0.
#' The null hypothesis is that the data is consistent with the DAG structure.
#'
#' @examples
#' # Create and test a chain DAG: X1 -> X2 -> X3
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' # Generate data consistent with the DAG
#' set.seed(123)
#' n <- 200
#' X1 <- rnorm(n)
#' X2 <- 0.5 * X1 + rnorm(n, sd = 0.5)
#' X3 <- 0.5 * X2 + rnorm(n, sd = 0.5)
#' data <- cbind(X1, X2, X3)
#'
#' # Test DAG compatibility
#' result <- mi_test(data, g, B = 100)
#' print(result)
#'
#' @export
mi_test <- function(data, dag,
                     test_functions = "squared",
                     ordering = "optimal",
                     B = 500,
                     alpha = 0.05,
                     parallel = FALSE,
                     verbose = TRUE,
                     seed = NULL) {
  # Validate inputs
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  n_obs <- nrow(data)
  n_vars <- ncol(data)

  if (!is.null(seed)) set.seed(seed)

  # Handle ordering selection
  if (is.character(ordering) && length(ordering) == 1) {
    if (ordering == "all") {
      return(.mi_test_all_orderings(data, dag, test_functions,
                                     B, alpha, parallel, verbose, seed))
    }

    if (ordering %in% c("optimal", "first", "random")) {
      ord_method <- ordering
      ord_result <- select_ordering(data, dag, method = ordering)
      ordering <- ord_result$ordering
      if (verbose && ord_method == "optimal") {
        cat("Selected ordering:", paste(ordering, collapse = " -> "), "\n")
        cat("Q_n =", format(ord_result$qn, digits = 4), "\n")
      }
    }
  }

  # Validate ordering
  if (!is_valid_ordering(dag, ordering)) {
    stop("Invalid topological ordering provided")
  }

  # Parse test functions
  h_spec <- .parse_test_functions(test_functions)

  if (verbose) {
    cat("Testing DAG with", n_obs, "observations and", n_vars, "variables\n")
    cat("Ordering:", paste(ordering, collapse = " -> "), "\n")
    cat("Test function(s):", paste(names(h_spec), collapse = ", "), "\n")
  }

  # Compute observed statistic
  qn_result <- compute_qn(data, dag, ordering = ordering)
  qn <- qn_result$qn

  if (verbose) {
    cat("Q_n^G =", format(qn, digits = 4), "\n")
  }

  # Compute test statistic for each test function
  test_results <- lapply(h_spec, function(spec) {
    compute_test_stat(data, dag, h = spec$h, h_params = spec$params,
                       qn_result = qn_result, ordering = ordering)
  })

  # Use first test function for bootstrap
  main_result <- test_results[[1]]
  statistic <- main_result$statistic

  if (verbose) {
    cat("Observed T_h =", format(statistic, digits = 4), "\n")
  }

  # Bootstrap inference
  if (B > 0) {
    if (verbose) cat("\nRunning bootstrap...\n")

    boot_result <- constrained_bootstrap(
      data, dag, B = B,
      h = h_spec[[1]]$h, h_params = h_spec[[1]]$params,
      ordering = ordering,
      parallel = parallel,
      verbose = verbose
    )

    p_value <- boot_result$p_value
  } else {
    boot_result <- NULL
    p_value <- NA
  }

  # Decision
  decision <- if (is.na(p_value)) {
    "not_computed"
  } else if (p_value < alpha) {
    "reject"
  } else {
    "fail_to_reject"
  }

  if (verbose && !is.na(p_value)) {
    cat("\np-value:", format(p_value, digits = 4), "\n")
    cat("Decision:", decision, "H0 at alpha =", alpha, "\n")
  }

  structure(
    list(
      statistic = statistic,
      p_value = p_value,
      decision = decision,
      alpha = alpha,
      qn = qn,
      lhs = main_result$lhs,
      rhs = main_result$rhs,
      bootstrap = boot_result,
      test_results = test_results,
      ordering = ordering,
      dag = dag,
      n_obs = n_obs,
      n_vars = n_vars
    ),
    class = "mi_test_result"
  )
}

#' Parse Test Function Specifications
#' @keywords internal
.parse_test_functions <- function(test_functions) {
  if (is.function(test_functions)) {
    return(list(custom = list(h = test_functions, params = list())))
  }

  if (is.character(test_functions) && length(test_functions) == 1) {
    return(setNames(list(list(h = test_functions, params = list())),
                    test_functions))
  }

  if (is.list(test_functions)) {
    result <- lapply(test_functions, function(spec) {
      if (is.function(spec)) {
        list(h = spec, params = list())
      } else if (is.character(spec)) {
        list(h = spec, params = list())
      } else if (is.list(spec)) {
        list(h = spec$h %||% "squared",
             params = spec$h_params %||% spec$params %||% list())
      }
    })

    if (is.null(names(result))) {
      names(result) <- paste0("h", seq_along(result))
    }

    return(result)
  }

  stop("Invalid test_functions specification")
}

#' Test All Orderings
#' @keywords internal
.mi_test_all_orderings <- function(data, dag, test_functions,
                                     B, alpha, parallel, verbose, seed) {
  orderings <- get_topological_orders(dag)
  n_orderings <- length(orderings)

  if (verbose) {
    cat("Testing", n_orderings, "topological orderings\n")
  }

  # Test each ordering
  results <- vector("list", n_orderings)
  p_values <- numeric(n_orderings)

  for (i in seq_len(n_orderings)) {
    if (verbose) {
      cat("\nOrdering", i, "/", n_orderings, ":",
          paste(orderings[[i]], collapse = " -> "), "\n")
    }

    results[[i]] <- mi_test(
      data, dag,
      test_functions = test_functions,
      ordering = orderings[[i]],
      B = B,
      alpha = alpha,  # Will adjust below
      parallel = parallel,
      verbose = FALSE,
      seed = if (!is.null(seed)) seed + i else NULL
    )

    p_values[i] <- results[[i]]$p_value
  }

  # Apply Bonferroni correction
  bonf <- bonferroni_orderings(p_values, alpha)

  # Also compute Fisher combination
  fisher <- fisher_combine(p_values)

  if (verbose) {
    cat("\n========================================\n")
    cat("Combined Results (", n_orderings, " orderings)\n", sep = "")
    cat("----------------------------------------\n")
    cat("Minimum p-value:", format(bonf$min_p, digits = 4), "\n")
    cat("Bonferroni-adjusted:", format(bonf$adjusted_min_p, digits = 4), "\n")
    cat("Fisher combined p:", format(fisher$combined_p, digits = 4), "\n")
    cat("Reject H0 (Bonferroni):", bonf$global_reject, "\n")
  }

  structure(
    list(
      results = results,
      p_values = p_values,
      orderings = orderings,
      bonferroni = bonf,
      fisher = fisher,
      global_decision = if (bonf$global_reject) "reject" else "fail_to_reject",
      alpha = alpha,
      dag = dag,
      n_orderings = n_orderings
    ),
    class = c("mi_test_all_orderings", "mi_test_result")
  )
}

#' Print Method for mi_test_result
#' @export
print.mi_test_result <- function(x, ...) {
  cat("\n")
  cat("Multilinear Inequality Test for DAG Compatibility\n")
  cat("==================================================\n\n")

  cat("DAG:", x$dag$n, "nodes,",
      sum(x$dag$adjacency), "edges\n")
  cat("Sample size:", x$n_obs, "\n")
  cat("Ordering:", paste(x$ordering, collapse = " -> "), "\n\n")

  cat("Test Results:\n")
  cat("  Q_n^G:", format(x$qn, digits = 4), "\n")
  cat("  T_h^G:", format(x$statistic, digits = 4), "\n")
  cat("    LHS (observed):", format(x$lhs, digits = 4), "\n")
  cat("    RHS (bound):", format(x$rhs, digits = 4), "\n")

  if (!is.na(x$p_value)) {
    cat("  p-value:", format(x$p_value, digits = 4), "\n\n")
    cat("Decision: ", x$decision, " H0 at alpha = ", x$alpha, "\n", sep = "")

    interpretation <- if (x$decision == "reject") {
      "Evidence against DAG compatibility"
    } else {
      "No evidence against DAG compatibility"
    }
    cat("Interpretation:", interpretation, "\n")
  } else {
    cat("  p-value: not computed (B = 0)\n")
  }

  invisible(x)
}

#' Print Method for All-Orderings Test
#' @export
print.mi_test_all_orderings <- function(x, ...) {
  cat("\n")
  cat("Multilinear Inequality Test - All Orderings\n")
  cat("============================================\n\n")

  cat("DAG:", x$dag$n, "nodes\n")
  cat("Number of orderings tested:", x$n_orderings, "\n\n")

  cat("P-values by ordering:\n")
  for (i in seq_along(x$orderings)) {
    cat("  ", paste(x$orderings[[i]], collapse = "->"), ": ",
        format(x$p_values[i], digits = 4), "\n", sep = "")
  }

  cat("\nCombined inference:\n")
  cat("  Minimum p-value:", format(x$bonferroni$min_p, digits = 4), "\n")
  cat("  Bonferroni-adjusted:", format(x$bonferroni$adjusted_min_p, digits = 4), "\n")
  cat("  Fisher combined:", format(x$fisher$combined_p, digits = 4), "\n\n")

  cat("Decision:", x$global_decision, "H0 at alpha =", x$alpha, "\n")

  invisible(x)
}

#' Plot Method for mi_test_result
#' @export
plot.mi_test_result <- function(x, type = "bootstrap", ...) {
  if (!is.null(x$bootstrap)) {
    plot(x$bootstrap, ...)
  } else {
    message("No bootstrap results to plot (B = 0)")
  }
  invisible(x)
}

#' Compare Two DAGs
#'
#' @description
#' Tests which of two DAGs better fits the observed data.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag1 First dag object
#' @param dag2 Second dag object
#' @param B Number of bootstrap replicates
#' @param alpha Significance level
#' @param parallel Use parallel processing
#' @param verbose Show progress
#'
#' @return A list with test results for both DAGs and comparison
#'
#' @examples
#' # Chain vs fork
#' A1 <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' A2 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
#'
#' g1 <- dag(A1, nodes = c("X1", "X2", "X3"))
#' g2 <- dag(A2, nodes = c("X1", "X2", "X3"))
#'
#' set.seed(123)
#' X1 <- rnorm(200)
#' X2 <- 0.5 * X1 + rnorm(200, sd = 0.5)
#' X3 <- 0.5 * X2 + rnorm(200, sd = 0.5)
#' data <- cbind(X1, X2, X3)
#'
#' compare_dags(data, g1, g2, B = 100)
#'
#' @export
compare_dags <- function(data, dag1, dag2,
                          B = 500, alpha = 0.05,
                          parallel = FALSE, verbose = TRUE) {
  if (verbose) cat("Testing DAG 1...\n")
  result1 <- mi_test(data, dag1, B = B, alpha = alpha,
                      parallel = parallel, verbose = verbose)

  if (verbose) cat("\nTesting DAG 2...\n")
  result2 <- mi_test(data, dag2, B = B, alpha = alpha,
                      parallel = parallel, verbose = verbose)

  # Comparison
  comparison <- list(
    dag1_statistic = result1$statistic,
    dag2_statistic = result2$statistic,
    dag1_p = result1$p_value,
    dag2_p = result2$p_value,
    dag1_qn = result1$qn,
    dag2_qn = result2$qn
  )

  # Preference based on p-values
  if (!is.na(result1$p_value) && !is.na(result2$p_value)) {
    if (result1$p_value > result2$p_value) {
      comparison$preferred <- "dag1"
      comparison$reason <- "Higher p-value (less evidence against)"
    } else if (result2$p_value > result1$p_value) {
      comparison$preferred <- "dag2"
      comparison$reason <- "Higher p-value (less evidence against)"
    } else {
      comparison$preferred <- "neither"
      comparison$reason <- "Equal p-values"
    }
  }

  structure(
    list(
      result1 = result1,
      result2 = result2,
      comparison = comparison,
      dag1 = dag1,
      dag2 = dag2
    ),
    class = "dag_comparison"
  )
}

#' Print Method for DAG Comparison
#' @export
print.dag_comparison <- function(x, ...) {
  cat("\n")
  cat("DAG Comparison Results\n")
  cat("======================\n\n")

  cat("DAG 1:\n")
  cat("  T_h:", format(x$comparison$dag1_statistic, digits = 4), "\n")
  cat("  p-value:", format(x$comparison$dag1_p, digits = 4), "\n")
  cat("  Q_n:", format(x$comparison$dag1_qn, digits = 4), "\n\n")

  cat("DAG 2:\n")
  cat("  T_h:", format(x$comparison$dag2_statistic, digits = 4), "\n")
  cat("  p-value:", format(x$comparison$dag2_p, digits = 4), "\n")
  cat("  Q_n:", format(x$comparison$dag2_qn, digits = 4), "\n\n")

  cat("Comparison:\n")
  cat("  Preferred:", x$comparison$preferred, "\n")
  cat("  Reason:", x$comparison$reason, "\n")

  invisible(x)
}
