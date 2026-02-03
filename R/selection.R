#' Test Function Selection and Diagnostics
#'
#' @description
#' Adaptive test function selection and diagnostic utilities for the
#' multilinear inequality test.
#'
#' @name selection
NULL

#' Adaptive Test Function Selection
#'
#' @description
#' Selects the optimal test function by data splitting, as described in
#' Algorithm 4 of Zambrano (2026). The data is split into a selection set
#' (to choose the test function) and an inference set (to compute the
#' final p-value), ensuring valid inference.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object specifying the DAG structure
#' @param candidates Named list of test function specifications. Each element
#'   should be a list with components \code{h} (function name or function) and
#'   optionally \code{h_params} (list of parameters). If NULL, uses default
#'   candidates.
#' @param split_ratio Proportion of data used for selection (default 0.5)
#' @param ordering Topological ordering to use (default "optimal")
#' @param B Number of bootstrap replicates for final inference (default 500)
#' @param alpha Significance level (default 0.05)
#' @param verbose Logical; show progress (default TRUE)
#' @param seed Random seed for reproducibility
#'
#' @return An S3 object of class "select_test_result" with components:
#'   \item{selected}{Name of selected test function}
#'   \item{selection_stats}{Test statistics on selection set for each candidate}
#'   \item{final_result}{mi_test result using selected function on inference set}
#'   \item{candidates}{The candidate test functions evaluated}
#'   \item{split_ratio}{The split ratio used}
#'   \item{n_selection}{Number of observations in selection set}
#'   \item{n_inference}{Number of observations in inference set}
#'
#' @details
#' The adaptive selection procedure works as follows:
#' \enumerate{
#'   \item Split data into selection set D_s and inference set D_i
#'   \item For each candidate test function h, compute T_h on D_s
#'   \item Select h* = argmax_h T_h (larger T_h indicates more power)
#'   \item Compute final p-value using h* on D_i only
#' }
#'
#' This data-splitting approach ensures that the final p-value is valid
#' (not inflated by selection), while still benefiting from adaptive
#' test function choice.
#'
#' @examples
#' # Create chain DAG
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' # Generate data
#' set.seed(123)
#' n <- 400
#' X1 <- rnorm(n)
#' X2 <- 0.5 * X1 + rnorm(n, sd = 0.5)
#' X3 <- 0.5 * X2 + rnorm(n, sd = 0.5)
#' data <- cbind(X1, X2, X3)
#'
#' # Adaptive selection
#' result <- select_test_function(data, g, B = 100)
#' print(result)
#'
#' @export
select_test_function <- function(data, dag,
                                  candidates = NULL,
                                  split_ratio = 0.5,
                                  ordering = "optimal",
                                  B = 500,
                                  alpha = 0.05,
                                  verbose = TRUE,
                                  seed = NULL) {
  # Validate inputs
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  n_obs <- nrow(data)
  n_vars <- ncol(data)

  if (split_ratio <= 0 || split_ratio >= 1) {
    stop("'split_ratio' must be between 0 and 1")
  }

  if (!is.null(seed)) set.seed(seed)

  # Default candidate test functions
  if (is.null(candidates)) {
    candidates <- list(
      squared = list(h = "squared", h_params = list()),
      poly4 = list(h = "poly", h_params = list(k = 2)),
      poly6 = list(h = "poly", h_params = list(k = 3)),
      exp025 = list(h = "exp", h_params = list(t = 0.25)),
      exp05 = list(h = "exp", h_params = list(t = 0.5)),
      indicator0 = list(h = "indicator", h_params = list(c = 0))
    )
  }

  # Step 1: Split data
  n_selection <- floor(n_obs * split_ratio)
  n_inference <- n_obs - n_selection

  if (n_selection < 50 || n_inference < 50) {
    warning("Small split sizes may lead to unstable results. ",
            "Consider using more data or a different split ratio.")
  }

  idx_selection <- sample(n_obs, n_selection)
  idx_inference <- setdiff(seq_len(n_obs), idx_selection)

  data_selection <- data[idx_selection, , drop = FALSE]
  data_inference <- data[idx_inference, , drop = FALSE]

  if (verbose) {
    cat("Adaptive Test Function Selection\n")
    cat("=================================\n")
    cat("Total observations:", n_obs, "\n")
    cat("Selection set:", n_selection, "observations\n")
    cat("Inference set:", n_inference, "observations\n")
    cat("Candidates:", length(candidates), "test functions\n\n")
  }

  # Step 2: Compute ordering on selection data
  if (is.character(ordering) && length(ordering) == 1 &&
      ordering %in% c("optimal", "first", "random")) {
    ord_result <- select_ordering(data_selection, dag, method = ordering)
    ordering <- ord_result$ordering
    if (verbose) {
      cat("Selected ordering:", paste(ordering, collapse = " -> "), "\n\n")
    }
  }

  # Step 3: Compute Q_n on selection data
  qn_result_sel <- compute_qn(data_selection, dag, ordering = ordering)

  # Step 4: Evaluate each candidate on selection set
  selection_stats <- numeric(length(candidates))
  names(selection_stats) <- names(candidates)

  if (verbose) cat("Evaluating candidates on selection set:\n")

  for (i in seq_along(candidates)) {
    cand <- candidates[[i]]
    cand_name <- names(candidates)[i]

    # Check admissibility
    h_func <- .parse_test_function(cand$h, cand$h_params %||% list())
    admiss <- check_admissibility(h_func, data_selection)

    if (!admiss$admissible) {
      selection_stats[i] <- -Inf
      if (verbose) {
        cat("  ", cand_name, ": INADMISSIBLE (", admiss$message, ")\n", sep = "")
      }
      next
    }

    # Compute test statistic
    result <- compute_test_stat(
      data_selection, dag,
      h = cand$h,
      h_params = cand$h_params %||% list(),
      qn_result = qn_result_sel,
      ordering = ordering
    )

    selection_stats[i] <- result$statistic

    if (verbose) {
      cat("  ", cand_name, ": T_h = ", format(result$statistic, digits = 4), "\n", sep = "")
    }
  }

  # Step 5: Select best candidate (maximum T_h)
  # Larger T_h = RHS - LHS means more "room" before violation, suggesting better power
  valid_stats <- selection_stats[is.finite(selection_stats)]

  if (length(valid_stats) == 0) {
    stop("No admissible test functions found")
  }

  selected_name <- names(which.max(valid_stats))
  selected_spec <- candidates[[selected_name]]

  if (verbose) {
    cat("\nSelected test function:", selected_name, "\n\n")
  }

  # Step 6: Compute final p-value on inference set
  if (verbose) cat("Computing p-value on inference set...\n")

  # Recompute ordering for inference set to be conservative
  # (or use same ordering for consistency - using same for now)
  final_result <- mi_test(
    data_inference, dag,
    test_functions = list(
      list(h = selected_spec$h, h_params = selected_spec$h_params %||% list())
    ),
    ordering = ordering,
    B = B,
    alpha = alpha,
    verbose = verbose,
    seed = if (!is.null(seed)) seed + 1000 else NULL
  )

  structure(
    list(
      selected = selected_name,
      selected_spec = selected_spec,
      selection_stats = selection_stats,
      final_result = final_result,
      candidates = candidates,
      split_ratio = split_ratio,
      n_selection = n_selection,
      n_inference = n_inference,
      ordering = ordering
    ),
    class = "select_test_result"
  )
}

#' Print Method for select_test_result
#' @export
print.select_test_result <- function(x, ...) {
  cat("\n")
  cat("Adaptive Test Function Selection Results\n")
  cat("=========================================\n\n")

  cat("Data Split:\n")
  cat("  Selection set:", x$n_selection, "observations\n")
  cat("  Inference set:", x$n_inference, "observations\n\n")

  cat("Selection Statistics (on selection set):\n")
  for (i in seq_along(x$selection_stats)) {
    marker <- if (names(x$selection_stats)[i] == x$selected) " <-- SELECTED" else ""
    val <- if (is.finite(x$selection_stats[i])) {
      format(x$selection_stats[i], digits = 4)
    } else {
      "inadmissible"
    }
    cat("  ", names(x$selection_stats)[i], ": T_h = ", val, marker, "\n", sep = "")
  }

  cat("\nFinal Inference (on inference set):\n")
  cat("  Test function:", x$selected, "\n")
  cat("  T_h:", format(x$final_result$statistic, digits = 4), "\n")
  cat("  p-value:", format(x$final_result$p_value, digits = 4), "\n")
  cat("  Decision:", x$final_result$decision, "\n")

  invisible(x)
}


#' Power Computation for Multilinear Inequality Test
#'
#' @description
#' Computes statistical power of the MI test via simulation. This function
#' helps with sample size planning by estimating the probability of correctly
#' rejecting an incorrect DAG hypothesis.
#'
#' @param dag_true True DAG generating the data
#' @param dag_test DAG hypothesis to test (default: same as dag_true)
#' @param n Sample size
#' @param effect_size Coefficient magnitude for edges (default 0.5)
#' @param noise_sd Standard deviation of noise terms (default 1)
#' @param n_sim Number of simulation replicates (default 200)
#' @param B Bootstrap replicates per simulation (default 200)
#' @param alpha Significance level (default 0.05)
#' @param test_function Test function to use (default "squared")
#' @param h_params Parameters for test function
#' @param ordering Ordering method (default "optimal")
#' @param parallel Use parallel processing (default FALSE)
#' @param verbose Show progress (default TRUE)
#' @param seed Random seed for reproducibility
#'
#' @return An S3 object of class "power_result" with components:
#'   \item{power}{Estimated power (rejection rate)}
#'   \item{power_se}{Standard error of power estimate}
#'   \item{power_ci}{95\% confidence interval for power}
#'   \item{rejection_rate}{Same as power}
#'   \item{p_values}{Vector of p-values from each simulation}
#'   \item{statistics}{Vector of test statistics from each simulation}
#'   \item{n}{Sample size used}
#'   \item{n_sim}{Number of simulations}
#'   \item{alpha}{Significance level}
#'   \item{dag_true}{True DAG}
#'   \item{dag_test}{Tested DAG}
#'
#' @details
#' For power analysis, use \code{dag_true} != \code{dag_test} to estimate
#' power to detect DAG misspecification. For Type I error calibration,
#' use \code{dag_true} == \code{dag_test}.
#'
#' Data is generated from a linear structural equation model:
#' X_j = sum_{k in pa(j)} beta_{kj} * X_k + epsilon_j
#' where epsilon_j ~ N(0, noise_sd^2) and beta coefficients have magnitude
#' \code{effect_size}.
#'
#' @examples
#' # Chain vs Fork power comparison
#' A_chain <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' A_fork <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
#'
#' chain <- dag(A_chain, nodes = c("X1", "X2", "X3"))
#' fork <- dag(A_fork, nodes = c("X1", "X2", "X3"))
#'
#' # Power to detect fork when true model is chain
#' # (Use small n_sim for example)
#' power_result <- power_mi_test(
#'   dag_true = chain,
#'   dag_test = fork,
#'   n = 300,
#'   n_sim = 50,
#'   B = 100
#' )
#' print(power_result)
#'
#' @export
power_mi_test <- function(dag_true,
                           dag_test = dag_true,
                           n = 200,
                           effect_size = 0.5,
                           noise_sd = 1,
                           n_sim = 200,
                           B = 200,
                           alpha = 0.05,
                           test_function = "squared",
                           h_params = list(),
                           ordering = "optimal",
                           parallel = FALSE,
                           verbose = TRUE,
                           seed = NULL) {
  # Validate inputs
  if (!is_dag(dag_true)) stop("'dag_true' must be a dag object")
  if (!is_dag(dag_test)) stop("'dag_test' must be a dag object")

  if (dag_true$n != dag_test$n) {
    stop("dag_true and dag_test must have the same number of nodes")
  }

  if (n < 50) {
    warning("Small sample size (n < 50) may lead to unreliable results")
  }

  if (!is.null(seed)) set.seed(seed)

  n_vars <- dag_true$n
  nodes <- dag_true$nodes

  # Determine if this is Type I error or power analysis
  same_dag <- identical(dag_true$adjacency, dag_test$adjacency)

  if (verbose) {
    cat("Power Analysis for Multilinear Inequality Test\n")
    cat("===============================================\n")
    if (same_dag) {
      cat("Mode: Type I error calibration (H0 is true)\n")
    } else {
      cat("Mode: Power analysis (H0 is false)\n")
    }
    cat("Sample size: n =", n, "\n")
    cat("Simulations:", n_sim, "\n")
    cat("Bootstrap replicates:", B, "\n")
    cat("Effect size:", effect_size, "\n\n")
  }

  # Storage for results
  p_values <- numeric(n_sim)
  statistics <- numeric(n_sim)

  # Progress tracking
  if (verbose) {
    cat("Running simulations...\n")
    pb_interval <- max(1, floor(n_sim / 20))
  }

  for (sim in seq_len(n_sim)) {
    if (verbose && sim %% pb_interval == 0) {
      cat(sprintf("\r  Progress: %d%%", round(100 * sim / n_sim)))
    }

    # Generate data from dag_true
    data <- .generate_dag_data(dag_true, n, effect_size, noise_sd,
                                seed = if (!is.null(seed)) seed + sim else NULL)

    # Test dag_test
    result <- tryCatch({
      mi_test(data, dag_test,
              test_functions = list(list(h = test_function, h_params = h_params)),
              ordering = ordering,
              B = B,
              alpha = alpha,
              verbose = FALSE)
    }, error = function(e) {
      list(p_value = NA, statistic = NA)
    })

    p_values[sim] <- result$p_value
    statistics[sim] <- result$statistic
  }

  if (verbose) cat("\r  Progress: 100%\n\n")

  # Compute power/Type I error
  valid_p <- p_values[!is.na(p_values)]
  n_valid <- length(valid_p)

  if (n_valid < n_sim * 0.9) {
    warning(sprintf("%.0f%% of simulations failed", 100 * (1 - n_valid / n_sim)))
  }

  rejection_rate <- mean(valid_p < alpha)
  power_se <- sqrt(rejection_rate * (1 - rejection_rate) / n_valid)
  power_ci <- c(
    max(0, rejection_rate - 1.96 * power_se),
    min(1, rejection_rate + 1.96 * power_se)
  )

  if (verbose) {
    cat("Results:\n")
    if (same_dag) {
      cat("  Type I error rate:", format(rejection_rate, digits = 3), "\n")
      cat("  (Nominal alpha:", alpha, ")\n")
    } else {
      cat("  Estimated power:", format(rejection_rate, digits = 3), "\n")
    }
    cat("  95% CI: [", format(power_ci[1], digits = 3), ", ",
        format(power_ci[2], digits = 3), "]\n", sep = "")
    cat("  Valid simulations:", n_valid, "/", n_sim, "\n")
  }

  structure(
    list(
      power = rejection_rate,
      power_se = power_se,
      power_ci = power_ci,
      rejection_rate = rejection_rate,
      p_values = p_values,
      statistics = statistics,
      n = n,
      n_sim = n_sim,
      n_valid = n_valid,
      B = B,
      alpha = alpha,
      effect_size = effect_size,
      noise_sd = noise_sd,
      test_function = test_function,
      dag_true = dag_true,
      dag_test = dag_test,
      same_dag = same_dag
    ),
    class = "power_result"
  )
}

#' Generate Data from DAG
#' @keywords internal
.generate_dag_data <- function(dag, n, effect_size = 0.5, noise_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_vars <- dag$n
  nodes <- dag$nodes
  A <- dag$adjacency

  # Get topological order
  order <- topological_orders(dag)[[1]]

  # Generate data
  data <- matrix(NA, nrow = n, ncol = n_vars)
  colnames(data) <- nodes

  for (j in order) {
    j_idx <- which(nodes == j)
    pa_j <- parents(dag, j)

    if (length(pa_j) == 0) {
      # Root node
      data[, j_idx] <- rnorm(n, sd = noise_sd)
    } else {
      # Child node: linear combination of parents
      pa_idx <- match(pa_j, nodes)
      data[, j_idx] <- rnorm(n, sd = noise_sd)
      for (pa in pa_idx) {
        data[, j_idx] <- data[, j_idx] + effect_size * data[, pa]
      }
    }
  }

  data
}

#' Print Method for power_result
#' @export
print.power_result <- function(x, ...) {
  cat("\n")
  cat("Power Analysis Results\n")
  cat("======================\n\n")

  if (x$same_dag) {
    cat("Mode: Type I Error Calibration\n")
    cat("  Estimated Type I error:", format(x$power, digits = 3), "\n")
    cat("  Nominal alpha:", x$alpha, "\n")
    diff <- x$power - x$alpha
    if (abs(diff) > 2 * x$power_se) {
      if (diff > 0) {
        cat("  WARNING: Test appears anti-conservative\n")
      } else {
        cat("  Test appears conservative\n")
      }
    } else {
      cat("  Type I error is well-calibrated\n")
    }
  } else {
    cat("Mode: Power Analysis\n")
    cat("  Estimated power:", format(x$power, digits = 3), "\n")
  }

  cat("  95% CI: [", format(x$power_ci[1], digits = 3), ", ",
      format(x$power_ci[2], digits = 3), "]\n\n", sep = "")

  cat("Settings:\n")
  cat("  Sample size: n =", x$n, "\n")
  cat("  Effect size:", x$effect_size, "\n")
  cat("  Simulations:", x$n_valid, "/", x$n_sim, "valid\n")
  cat("  Bootstrap replicates:", x$B, "\n")
  cat("  Test function:", x$test_function, "\n")

  invisible(x)
}


#' Diagnose Potential Confounding or Misspecification
#'
#' @description
#' Provides diagnostic information when test results are unexpected,
#' particularly when all candidate DAGs are rejected. This may indicate
#' unmeasured confounding, model misspecification, or other issues.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag_list List of dag objects to test
#' @param B Number of bootstrap replicates (default 300)
#' @param alpha Significance level (default 0.05)
#' @param test_function Test function to use (default "squared")
#' @param verbose Show progress (default TRUE)
#'
#' @return An S3 object of class "confounding_diagnostic" with components:
#'   \item{all_rejected}{Logical; TRUE if all DAGs were rejected}
#'   \item{results}{List of mi_test results for each DAG}
#'   \item{p_values}{Vector of p-values}
#'   \item{min_p}{Minimum p-value across DAGs}
#'   \item{best_dag}{Index of DAG with highest p-value}
#'   \item{diagnosis}{Character string describing likely issue}
#'   \item{recommendations}{Character vector of recommended actions}
#'   \item{residual_analysis}{Optional residual diagnostics}
#'
#' @details
#' The diagnostic procedure:
#' \enumerate{
#'   \item Tests all candidate DAGs
#'   \item If ALL are rejected, flags potential confounding
#'   \item Analyzes residual patterns from best-fitting DAG
#'   \item Provides actionable recommendations
#' }
#'
#' Common causes of all-DAG rejection:
#' \itemize{
#'   \item Unmeasured confounding (latent common causes)
#'   \item Selection bias in the sample
#'   \item Measurement error
#'   \item Non-causal associations (e.g., reverse causation)
#'   \item Model misspecification (e.g., nonlinear relationships)
#' }
#'
#' @examples
#' # Create candidate DAGs for 3-variable system
#' A1 <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' A2 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
#' A3 <- matrix(c(0, 0, 0, 1, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#'
#' chain <- dag(A1, nodes = c("X1", "X2", "X3"))
#' fork <- dag(A2, nodes = c("X1", "X2", "X3"))
#' collider <- dag(A3, nodes = c("X1", "X2", "X3"))
#'
#' # Generate data with unmeasured confounder
#' set.seed(123)
#' n <- 300
#' U <- rnorm(n)  # Unmeasured confounder
#' X1 <- 0.5 * U + rnorm(n)
#' X2 <- 0.5 * U + rnorm(n)
#' X3 <- 0.5 * X2 + rnorm(n)
#' data <- cbind(X1, X2, X3)
#'
#' # Diagnose
#' diag <- diagnose_confounding(data, list(chain, fork, collider), B = 100)
#' print(diag)
#'
#' @export
diagnose_confounding <- function(data, dag_list,
                                   B = 300,
                                   alpha = 0.05,
                                   test_function = "squared",
                                   verbose = TRUE) {
  if (!is.list(dag_list)) {
    dag_list <- list(dag_list)
  }

  n_dags <- length(dag_list)

  if (n_dags == 0) {
    stop("dag_list must contain at least one DAG")
  }

  # Validate all DAGs have same nodes
  nodes <- dag_list[[1]]$nodes
  n_vars <- dag_list[[1]]$n
  for (i in seq_along(dag_list)) {
    if (!is_dag(dag_list[[i]])) {
      stop("Element ", i, " of dag_list is not a dag object")
    }
    if (!identical(dag_list[[i]]$nodes, nodes)) {
      stop("All DAGs must have the same nodes")
    }
  }

  data <- validate_data(data, dag_list[[1]])
  n_obs <- nrow(data)

  if (verbose) {
    cat("Confounding Diagnostic Analysis\n")
    cat("================================\n")
    cat("Testing", n_dags, "candidate DAG(s)\n")
    cat("Sample size: n =", n_obs, "\n\n")
  }

  # Test each DAG
  results <- vector("list", n_dags)
  p_values <- numeric(n_dags)

  for (i in seq_len(n_dags)) {
    if (verbose) {
      cat("Testing DAG", i, "/", n_dags, "...\n")
    }

    results[[i]] <- mi_test(
      data, dag_list[[i]],
      test_functions = test_function,
      B = B,
      alpha = alpha,
      verbose = FALSE
    )

    p_values[i] <- results[[i]]$p_value

    if (verbose) {
      cat("  p-value:", format(p_values[i], digits = 4), "\n")
    }
  }

  # Analysis
  all_rejected <- all(p_values < alpha)
  none_rejected <- all(p_values >= alpha)
  best_dag_idx <- which.max(p_values)
  min_p <- min(p_values)
  max_p <- max(p_values)

  # Diagnosis
  if (all_rejected) {
    diagnosis <- "ALL_REJECTED"
    severity <- "high"
    description <- paste0(
      "All ", n_dags, " candidate DAG(s) were rejected at alpha = ", alpha, ". ",
      "This strongly suggests model misspecification or unmeasured confounding."
    )
  } else if (none_rejected) {
    diagnosis <- "NONE_REJECTED"
    severity <- "low"
    description <- paste0(
      "No candidate DAGs were rejected. ",
      "The data is consistent with multiple causal structures."
    )
  } else {
    n_rejected <- sum(p_values < alpha)
    diagnosis <- "SOME_REJECTED"
    severity <- "medium"
    description <- paste0(
      n_rejected, " of ", n_dags, " DAG(s) rejected. ",
      "Consider focusing on non-rejected structure(s)."
    )
  }

  # Recommendations
  recommendations <- character(0)

  if (all_rejected) {
    recommendations <- c(
      "1. Consider unmeasured confounding: Are there latent common causes?",
      "2. Try FCI algorithm to search for latent variable structures",
      "3. Check for selection bias in data collection",
      "4. Examine potential measurement error",
      "5. Consider nonlinear relationships (try different test functions)",
      "6. Perform sensitivity analysis with different preprocessing"
    )

    # Additional check: examine correlations
    cors <- cor(data)
    high_cors <- which(abs(cors) > 0.7 & upper.tri(cors), arr.ind = TRUE)
    if (nrow(high_cors) > 0) {
      recommendations <- c(
        recommendations,
        paste0("7. Note: High correlations detected between variables ",
               "(may indicate strong confounding or near-collinearity)")
      )
    }
  } else if (none_rejected) {
    recommendations <- c(
      "1. All DAGs are consistent with data - cannot discriminate",
      "2. Consider collecting more data to increase power",
      "3. Try different test functions for better discrimination",
      "4. Focus on theoretical or domain knowledge to select among structures"
    )
  } else {
    recommendations <- c(
      "1. Focus analysis on non-rejected DAG structure(s)",
      "2. Report p-values for transparency about uncertainty",
      "3. Consider sensitivity analysis with rejected DAGs"
    )
  }

  # Residual analysis for best DAG (basic)
  best_dag <- dag_list[[best_dag_idx]]
  residual_info <- NULL

  if (all_rejected && n_vars <= 5) {
    # Simple residual check
    residual_info <- .compute_residual_diagnostics(data, best_dag)
  }

  if (verbose) {
    cat("\n")
    cat("Diagnosis:", diagnosis, "\n")
    cat("Description:", description, "\n")
    cat("\nRecommendations:\n")
    for (rec in recommendations) {
      cat("  ", rec, "\n", sep = "")
    }
  }

  structure(
    list(
      all_rejected = all_rejected,
      none_rejected = none_rejected,
      results = results,
      p_values = p_values,
      min_p = min_p,
      max_p = max_p,
      best_dag = best_dag_idx,
      n_dags = n_dags,
      alpha = alpha,
      diagnosis = diagnosis,
      severity = severity,
      description = description,
      recommendations = recommendations,
      residual_analysis = residual_info,
      dag_list = dag_list
    ),
    class = "confounding_diagnostic"
  )
}

#' Compute Basic Residual Diagnostics
#' @keywords internal
.compute_residual_diagnostics <- function(data, dag) {
  n_obs <- nrow(data)
  n_vars <- ncol(data)
  nodes <- dag$nodes

  residuals <- matrix(NA, nrow = n_obs, ncol = n_vars)
  colnames(residuals) <- nodes
  r_squared <- numeric(n_vars)
  names(r_squared) <- nodes

  for (j in seq_len(n_vars)) {
    pa_j <- parents(dag, nodes[j])

    if (length(pa_j) == 0) {
      # Root node: residual is the variable itself (centered)
      residuals[, j] <- data[, j] - mean(data[, j])
      r_squared[j] <- 0
    } else {
      # Regress on parents
      pa_idx <- match(pa_j, nodes)
      X <- data[, pa_idx, drop = FALSE]
      y <- data[, j]

      fit <- lm(y ~ X)
      residuals[, j] <- residuals(fit)
      r_squared[j] <- summary(fit)$r.squared
    }
  }

  # Check residual normality
  normality_p <- apply(residuals, 2, function(r) {
    if (n_obs > 5000) {
      # Use subset for large samples
      r <- sample(r, 5000)
    }
    shapiro.test(r)$p.value
  })

  # Check residual correlations (should be independent if DAG is correct)
  resid_cors <- cor(residuals)

  list(
    residuals = residuals,
    r_squared = r_squared,
    normality_p = normality_p,
    residual_correlations = resid_cors
  )
}

#' Print Method for confounding_diagnostic
#' @export
print.confounding_diagnostic <- function(x, ...) {
  cat("\n")
  cat("Confounding Diagnostic Results\n")
  cat("==============================\n\n")

  cat("DAGs Tested:", x$n_dags, "\n")
  cat("P-values:\n")
  for (i in seq_len(x$n_dags)) {
    status <- if (x$p_values[i] < x$alpha) "REJECTED" else "not rejected"
    best <- if (i == x$best_dag) " (best)" else ""
    cat("  DAG ", i, ": p = ", format(x$p_values[i], digits = 4),
        " - ", status, best, "\n", sep = "")
  }

  cat("\n")
  cat("Diagnosis:", x$diagnosis, "\n")
  cat("Severity:", x$severity, "\n")
  cat("\nDescription:\n  ", x$description, "\n", sep = "")

  cat("\nRecommendations:\n")
  for (rec in x$recommendations) {
    cat("  ", rec, "\n", sep = "")
  }

  if (!is.null(x$residual_analysis)) {
    cat("\nResidual Analysis (best DAG):\n")
    cat("  R-squared by variable:\n")
    for (i in seq_along(x$residual_analysis$r_squared)) {
      cat("    ", names(x$residual_analysis$r_squared)[i], ": ",
          format(x$residual_analysis$r_squared[i], digits = 3), "\n", sep = "")
    }
  }

  invisible(x)
}
