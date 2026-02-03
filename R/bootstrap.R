#' Constrained Bootstrap Inference
#'
#' @description
#' Implements the constrained bootstrap procedure for inference under H_0
#' (data is consistent with DAG structure).
#'
#' @name bootstrap
NULL

#' Constrained Bootstrap for DAG Testing
#'
#' @description
#' Generates bootstrap samples under the null hypothesis that the data
#' is consistent with the specified DAG structure, and computes the
#' bootstrap distribution of the test statistic.
#'
#' @param data Numeric matrix with n observations and p variables
#' @param dag A dag object specifying the DAG structure
#' @param B Number of bootstrap replicates (default 500)
#' @param h Test function (default "squared")
#' @param h_params Parameters for test function
#' @param ordering Optional topological ordering
#' @param parallel Logical or number of cores for parallel processing
#' @param verbose Logical; show progress (default TRUE)
#' @param seed Random seed for reproducibility
#'
#' @return An object of class "bootstrap_result" with components:
#'   \item{t0}{Observed test statistic}
#'   \item{t_boot}{Vector of bootstrap test statistics}
#'   \item{p_value}{Bootstrap p-value}
#'   \item{B}{Number of bootstrap replicates}
#'   \item{qn}{Q_n^G value}
#'   \item{se_boot}{Bootstrap standard error}
#'   \item{ci_lower}{Lower confidence bound (2.5%)}
#'   \item{ci_upper}{Upper confidence bound (97.5%)}
#'
#' @details
#' The constrained bootstrap generates samples that respect the DAG
#' structure under H_0. For each bootstrap replicate:
#' 1. Sample X_1 from the marginal distribution
#' 2. For j = 2, ..., n: sample X_j from p(X_j | pa(X_j))
#'
#' The p-value is the proportion of bootstrap statistics at least as
#' extreme as the observed statistic.
#'
#' @examples
#' A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
#' g <- dag(A, nodes = c("X1", "X2", "X3"))
#'
#' set.seed(123)
#' X1 <- rnorm(100)
#' X2 <- 0.5 * X1 + rnorm(100, sd = 0.5)
#' X3 <- 0.5 * X2 + rnorm(100, sd = 0.5)
#' data <- cbind(X1, X2, X3)
#'
#' # Bootstrap test
#' boot_result <- constrained_bootstrap(data, g, B = 100)
#'
#' @export
constrained_bootstrap <- function(data, dag, B = 500,
                                    h = "squared", h_params = list(),
                                    ordering = NULL,
                                    parallel = FALSE,
                                    verbose = TRUE,
                                    seed = NULL) {
  # Validate inputs
  if (!is_dag(dag)) stop("'dag' must be a dag object")
  data <- validate_data(data, dag)

  n_obs <- nrow(data)
  n_vars <- ncol(data)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Compute observed statistic and Q_n
  qn_result <- compute_qn(data, dag, ordering = ordering)
  t0_result <- compute_test_stat(data, dag, h = h, h_params = h_params,
                                   qn_result = qn_result)
  t0 <- t0_result$statistic
  ordering <- qn_result$ordering

  # Reorder data
  data_ordered <- data[, ordering, drop = FALSE]

  # Fit conditional distributions for bootstrap
  cond_kdes <- .fit_conditional_dists(data_ordered, dag, ordering)

  # Run bootstrap
  if (parallel && check_package("future", "parallel bootstrap")) {
    t_boot <- .bootstrap_parallel(
      data_ordered, dag, ordering, cond_kdes, qn_result,
      h, h_params, B, n_obs, n_vars, verbose
    )
  } else {
    t_boot <- .bootstrap_sequential(
      data_ordered, dag, ordering, cond_kdes, qn_result,
      h, h_params, B, n_obs, n_vars, verbose
    )
  }

  # Compute p-value (one-sided: test for T_h too small)
  # Under H0, T_h should be >= 0, so we test if observed is "too negative"
  p_value <- mean(t_boot <= t0)

  # Bootstrap statistics
  se_boot <- sd(t_boot)
  ci <- quantile(t_boot, c(0.025, 0.975))

  structure(
    list(
      t0 = t0,
      t_boot = t_boot,
      p_value = p_value,
      B = B,
      qn = qn_result$qn,
      se_boot = se_boot,
      ci_lower = ci[1],
      ci_upper = ci[2],
      ordering = ordering,
      h_name = t0_result$h_name,
      n_obs = n_obs,
      n_vars = n_vars
    ),
    class = "bootstrap_result"
  )
}

#' Fit Conditional Distributions for Bootstrap
#' @keywords internal
.fit_conditional_dists <- function(data, dag, ordering) {
  n_vars <- ncol(data)

  # For each variable, fit conditional distribution given parents
  cond_kdes <- vector("list", n_vars)

  for (j in seq_len(n_vars)) {
    node_name <- ordering[j]
    parent_names <- parents(dag, node_name)

    if (length(parent_names) == 0) {
      # Marginal distribution
      cond_kdes[[j]] <- list(
        type = "marginal",
        kde = kde_marginal(data[, j])
      )
    } else {
      # Conditional distribution
      # Find parent indices in ordered data
      parent_idx <- match(parent_names, ordering)

      if (length(parent_idx) == 1) {
        # Single parent - use conditional KDE
        cond_kdes[[j]] <- list(
          type = "conditional",
          kde = kde_conditional(data[, j], data[, parent_idx]),
          parent_idx = parent_idx
        )
      } else {
        # Multiple parents - use residual approach
        # Fit linear model and use residual distribution
        y <- data[, j]
        X <- data[, parent_idx, drop = FALSE]

        fit <- lm(y ~ X)
        residuals <- residuals(fit)
        coeffs <- coef(fit)

        cond_kdes[[j]] <- list(
          type = "linear_conditional",
          coeffs = coeffs,
          residual_kde = kde_marginal(residuals),
          parent_idx = parent_idx
        )
      }
    }
  }

  cond_kdes
}

#' Generate One Bootstrap Sample
#' @keywords internal
.generate_bootstrap_sample <- function(n_obs, n_vars, cond_kdes) {
  boot_data <- matrix(NA, nrow = n_obs, ncol = n_vars)

  for (j in seq_len(n_vars)) {
    cond <- cond_kdes[[j]]

    if (cond$type == "marginal") {
      # Sample from marginal
      boot_data[, j] <- .sample_from_kde(cond$kde, n_obs)

    } else if (cond$type == "conditional") {
      # Sample from conditional given parent
      parent_values <- boot_data[, cond$parent_idx]
      boot_data[, j] <- .sample_conditional_kde(
        cond$kde, parent_values, n_obs
      )

    } else if (cond$type == "linear_conditional") {
      # Linear model + residual
      parent_values <- boot_data[, cond$parent_idx, drop = FALSE]
      X <- cbind(1, parent_values)
      mean_values <- X %*% cond$coeffs
      residuals <- .sample_from_kde(cond$residual_kde, n_obs)
      boot_data[, j] <- mean_values + residuals
    }
  }

  boot_data
}

#' Sample from Marginal KDE
#' @keywords internal
.sample_from_kde <- function(kde, n) {
  # Mixture of Gaussians sampling
  # Sample from data with added noise
  idx <- sample(length(kde$data), n, replace = TRUE)
  kde$data[idx] + rnorm(n, sd = kde$bandwidth)
}

#' Sample from Conditional KDE
#' @keywords internal
.sample_conditional_kde <- function(kde, x_cond, n) {
  samples <- numeric(n)

  for (i in seq_len(n)) {
    # Use rejection sampling
    samples[i] <- kde_cond_sample(kde, x_cond[i], n_samples = 1)
  }

  samples
}

#' Sequential Bootstrap
#' @keywords internal
.bootstrap_sequential <- function(data, dag, ordering, cond_kdes, qn_result,
                                    h, h_params, B, n_obs, n_vars, verbose) {
  t_boot <- numeric(B)

  if (verbose) {
    cat("Running", B, "bootstrap replicates...\n")
  }

  for (b in seq_len(B)) {
    # Generate bootstrap sample
    boot_data <- .generate_bootstrap_sample(n_obs, n_vars, cond_kdes)

    # Compute test statistic
    boot_result <- compute_test_stat(
      boot_data, dag, h = h, h_params = h_params,
      ordering = ordering
    )
    t_boot[b] <- boot_result$statistic

    if (verbose && b %% 50 == 0) {
      cat(sprintf("  Completed %d/%d\n", b, B))
    }
  }

  if (verbose) cat("Done.\n")

  t_boot
}

#' Parallel Bootstrap
#' @keywords internal
.bootstrap_parallel <- function(data, dag, ordering, cond_kdes, qn_result,
                                  h, h_params, B, n_obs, n_vars, verbose) {
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("furrr", quietly = TRUE)) {
    warning("Packages 'future' and 'furrr' needed for parallel; using sequential")
    return(.bootstrap_sequential(data, dag, ordering, cond_kdes, qn_result,
                                  h, h_params, B, n_obs, n_vars, verbose))
  }

  old_plan <- setup_parallel(TRUE)
  on.exit(restore_parallel(old_plan), add = TRUE)

  if (verbose) {
    cat("Running", B, "bootstrap replicates in parallel...\n")
  }

  # Function for each bootstrap replicate
  boot_func <- function(b) {
    boot_data <- .generate_bootstrap_sample(n_obs, n_vars, cond_kdes)
    boot_result <- compute_test_stat(
      boot_data, dag, h = h, h_params = h_params,
      ordering = ordering
    )
    boot_result$statistic
  }

  t_boot <- furrr::future_map_dbl(
    seq_len(B),
    boot_func,
    .options = furrr::furrr_options(seed = TRUE)
  )

  if (verbose) cat("Done.\n")

  t_boot
}

#' Print Method for bootstrap_result
#' @export
print.bootstrap_result <- function(x, ...) {
  cat("Constrained Bootstrap Results\n")
  cat("=============================\n")
  cat("Observed T_h:", format(x$t0, digits = 4), "\n")
  cat("Bootstrap p-value:", format(x$p_value, digits = 4), "\n")
  cat("Bootstrap SE:", format(x$se_boot, digits = 4), "\n")
  cat("95% CI: [", format(x$ci_lower, digits = 4), ", ",
      format(x$ci_upper, digits = 4), "]\n", sep = "")
  cat("B =", x$B, "replicates\n")
  cat("Test function:", x$h_name, "\n")
  invisible(x)
}

#' Plot Method for bootstrap_result
#' @export
plot.bootstrap_result <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .plot_bootstrap_base(x, ...)
    return(invisible(x))
  }

  df <- data.frame(t = x$t_boot)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 30, fill = "lightblue", color = "black") +
    ggplot2::geom_vline(xintercept = x$t0, color = "red", linewidth = 1.2,
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "gray40", linewidth = 0.8) +
    ggplot2::labs(
      title = "Bootstrap Distribution of T_h^G",
      subtitle = sprintf("Observed = %.4f, p-value = %.4f", x$t0, x$p_value),
      x = "Test Statistic",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(x)
}

#' Base R Plot for Bootstrap
#' @keywords internal
.plot_bootstrap_base <- function(x, ...) {
  hist(x$t_boot, breaks = 30, freq = FALSE,
       main = "Bootstrap Distribution of T_h^G",
       xlab = "Test Statistic", col = "lightblue")
  abline(v = x$t0, col = "red", lwd = 2, lty = 2)
  abline(v = 0, col = "gray40")
  legend("topright",
         legend = c(paste("Observed =", round(x$t0, 4)),
                    paste("p =", round(x$p_value, 4))),
         lty = c(2, NA), col = c("red", NA))
}

#' Calibrate Bootstrap p-value
#'
#' @description
#' Computes calibrated p-values using the double bootstrap or other
#' corrections for better finite-sample performance.
#'
#' @param boot_result A bootstrap_result object
#' @param data Original data matrix
#' @param dag A dag object
#' @param B2 Number of second-level bootstrap replicates
#' @param method Calibration method: "percentile" (default), "bca", or "double"
#'
#' @return List with calibrated p-value and confidence interval
#'
#' @export
calibrate_pvalue <- function(boot_result, data, dag, B2 = 100,
                              method = "percentile") {
  method <- match.arg(method, c("percentile", "bca", "double"))

  if (method == "percentile") {
    # Standard percentile method (already computed)
    return(list(
      p_value = boot_result$p_value,
      ci = c(boot_result$ci_lower, boot_result$ci_upper),
      method = method
    ))
  }

  if (method == "bca") {
    # BCa correction
    return(.bca_correction(boot_result, data, dag))
  }

  if (method == "double") {
    # Double bootstrap
    return(.double_bootstrap(boot_result, data, dag, B2))
  }
}

#' BCa Correction
#' @keywords internal
.bca_correction <- function(boot_result, data, dag) {
  t0 <- boot_result$t0
  t_boot <- boot_result$t_boot
  B <- length(t_boot)

  # Bias correction factor z0
  z0 <- qnorm(mean(t_boot < t0))

  # Acceleration factor a (using jackknife)
  n <- boot_result$n_obs
  jack_vals <- numeric(n)

  for (i in seq_len(n)) {
    jack_data <- data[-i, , drop = FALSE]
    jack_result <- compute_test_stat(jack_data, dag,
                                       h = boot_result$h_name,
                                       ordering = boot_result$ordering)
    jack_vals[i] <- jack_result$statistic
  }

  jack_mean <- mean(jack_vals)
  a <- sum((jack_mean - jack_vals)^3) /
    (6 * (sum((jack_mean - jack_vals)^2))^1.5)

  # Adjusted percentiles
  alpha <- c(0.025, 0.975)
  z_alpha <- qnorm(alpha)
  adj_alpha <- pnorm(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))

  ci <- quantile(t_boot, adj_alpha)

  # Adjusted p-value
  adj_p <- pnorm(z0 + (z0 + qnorm(boot_result$p_value)) /
                   (1 - a * (z0 + qnorm(boot_result$p_value))))

  list(
    p_value = adj_p,
    ci = ci,
    z0 = z0,
    a = a,
    method = "bca"
  )
}

#' Double Bootstrap
#' @keywords internal
.double_bootstrap <- function(boot_result, data, dag, B2) {
  warning("Double bootstrap not fully implemented; returning uncalibrated")

  list(
    p_value = boot_result$p_value,
    ci = c(boot_result$ci_lower, boot_result$ci_upper),
    method = "double (uncalibrated)"
  )
}
