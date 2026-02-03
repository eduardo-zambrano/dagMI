#' Visualization Functions
#'
#' @description
#' Plotting and visualization functions for DAG testing results.
#'
#' @name visualization
NULL

#' Plot Bootstrap Distribution
#'
#' @description
#' Plots the bootstrap distribution of the test statistic with the
#' observed value highlighted.
#'
#' @param boot_result A bootstrap_result or mi_test_result object
#' @param show_ci Logical; show confidence interval (default TRUE)
#' @param show_zero Logical; show vertical line at zero (default TRUE)
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object
#'
#' @export
plot_bootstrap <- function(boot_result, show_ci = TRUE, show_zero = TRUE, ...) {
  # Extract bootstrap results
  if (inherits(boot_result, "mi_test_result")) {
    if (is.null(boot_result$bootstrap)) {
      stop("No bootstrap results available")
    }
    boot_result <- boot_result$bootstrap
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot(boot_result)
    return(invisible(NULL))
  }

  df <- data.frame(t = boot_result$t_boot)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 30, fill = "steelblue", color = "white", alpha = 0.7
    ) +
    ggplot2::geom_density(color = "darkblue", linewidth = 0.8)

  # Add observed value
  p <- p + ggplot2::geom_vline(
    xintercept = boot_result$t0,
    color = "red", linewidth = 1.2, linetype = "dashed"
  )

  # Add zero line
  if (show_zero) {
    p <- p + ggplot2::geom_vline(
      xintercept = 0, color = "gray40", linewidth = 0.8
    )
  }

  # Add CI
  if (show_ci) {
    p <- p + ggplot2::geom_vline(
      xintercept = c(boot_result$ci_lower, boot_result$ci_upper),
      color = "darkgreen", linewidth = 0.6, linetype = "dotted"
    )
  }

  # Labels
  p <- p + ggplot2::labs(
    title = "Bootstrap Distribution of Test Statistic",
    subtitle = sprintf("T_h = %.4f, p = %.4f",
                        boot_result$t0, boot_result$p_value),
    x = expression(T[h]^G),
    y = "Density"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  p
}

#' Plot Q_n Values Across Orderings
#'
#' @description
#' Visualizes Q_n values for different topological orderings.
#'
#' @param qn_all_result Result from compute_qn_all_orderings
#' @param highlight_best Logical; highlight the best ordering (default TRUE)
#'
#' @return A ggplot object
#'
#' @export
plot_qn_orderings <- function(qn_all_result, highlight_best = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    barplot(qn_all_result$qn_values,
            main = "Q_n by Topological Ordering",
            las = 2, col = "steelblue")
    return(invisible(NULL))
  }

  df <- data.frame(
    ordering = names(qn_all_result$qn_values),
    qn = qn_all_result$qn_values,
    stringsAsFactors = FALSE
  )

  # Identify best
  df$is_best <- df$qn == max(df$qn)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = reorder(.data$ordering, .data$qn),
    y = .data$qn,
    fill = .data$is_best
  )) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("steelblue", "coral")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = expression(Q[n]^G ~ "by Topological Ordering"),
      x = "Ordering",
      y = expression(Q[n]^G)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  p
}

#' Power Curve Plot
#'
#' @description
#' Plots empirical power as a function of sample size or effect size.
#'
#' @param power_results Data frame with columns: x (e.g., sample size),
#'   power, and optionally method
#' @param x_var Name of x variable
#' @param x_label Label for x-axis
#' @param title Plot title
#'
#' @return A ggplot object
#'
#' @export
plot_power_curve <- function(power_results, x_var = "n", x_label = "Sample Size",
                              title = "Empirical Power") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot(power_results[[x_var]], power_results$power,
         type = "b", xlab = x_label, ylab = "Power", main = title)
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(power_results, ggplot2::aes(
    x = .data[[x_var]],
    y = .data$power
  ))

  # Add method grouping if present
  if ("method" %in% names(power_results)) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(color = .data$method, linetype = .data$method),
                          linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(color = .data$method, shape = .data$method),
                           size = 3)
  } else {
    p <- p +
      ggplot2::geom_line(color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(color = "steelblue", size = 3)
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Power"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  p
}

#' Plot DAG Comparison
#'
#' @description
#' Side-by-side visualization of two DAGs with test results.
#'
#' @param comparison A dag_comparison object from compare_dags
#'
#' @return A combined plot
#'
#' @export
plot_dag_comparison <- function(comparison) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    par(mfrow = c(1, 2))
    plot(comparison$dag1, main = paste("DAG 1: p =",
                                         round(comparison$comparison$dag1_p, 4)))
    plot(comparison$dag2, main = paste("DAG 2: p =",
                                         round(comparison$comparison$dag2_p, 4)))
    par(mfrow = c(1, 1))
    return(invisible(NULL))
  }

  # Create summary data frame
  df <- data.frame(
    dag = c("DAG 1", "DAG 2"),
    p_value = c(comparison$comparison$dag1_p,
                 comparison$comparison$dag2_p),
    qn = c(comparison$comparison$dag1_qn,
            comparison$comparison$dag2_qn),
    statistic = c(comparison$comparison$dag1_statistic,
                   comparison$comparison$dag2_statistic)
  )

  df$preferred <- df$dag == paste0("DAG ", ifelse(
    comparison$comparison$preferred == "dag1", "1",
    ifelse(comparison$comparison$preferred == "dag2", "2", "neither")
  ))

  # P-value comparison
  p1 <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$dag, y = .data$p_value, fill = .data$preferred
  )) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("gray70", "coral")) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(title = "P-values", y = "P-value", x = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Q_n comparison
  p2 <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$dag, y = .data$qn, fill = .data$preferred
  )) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("gray70", "coral")) +
    ggplot2::labs(title = expression(Q[n]^G), y = expression(Q[n]^G), x = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p1 + p2 + patchwork::plot_annotation(
      title = "DAG Comparison",
      subtitle = paste("Preferred:", comparison$comparison$preferred)
    )
  } else {
    print(p1)
    print(p2)
  }
}

#' Subgraph Test Heatmap
#'
#' @description
#' Visualizes p-values across subgraphs as a heatmap.
#'
#' @param subgraph_result A subgraph_test_result object
#' @param show_adjusted Logical; use adjusted p-values (default TRUE)
#'
#' @return A ggplot object
#'
#' @export
plot_subgraph_pvalues <- function(subgraph_result, show_adjusted = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    p_vals <- if (show_adjusted) subgraph_result$adjusted_p else subgraph_result$p_values
    barplot(p_vals, names.arg = seq_along(p_vals),
            main = "Subgraph P-values", ylab = "P-value")
    abline(h = subgraph_result$alpha, lty = 2)
    return(invisible(NULL))
  }

  df <- data.frame(
    subgraph = sapply(subgraph_result$subgraph_nodes, paste, collapse = ","),
    p_value = if (show_adjusted) subgraph_result$adjusted_p else subgraph_result$p_values,
    rejected = subgraph_result$rejected,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = reorder(.data$subgraph, -.data$p_value),
    y = .data$p_value,
    fill = .data$rejected
  )) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = subgraph_result$alpha,
                         linetype = "dashed", color = "red") +
    ggplot2::scale_fill_manual(
      values = c("steelblue", "coral"),
      labels = c("Not Rejected", "Rejected"),
      name = "Decision"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Subgraph Test Results",
      subtitle = paste("Correction:", subgraph_result$correction),
      x = "Subgraph",
      y = ifelse(show_adjusted, "Adjusted P-value", "P-value")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  p
}

#' Diagnostic Plot for Test Assumptions
#'
#' @description
#' Creates diagnostic plots to check test assumptions.
#'
#' @param data Numeric matrix
#' @param dag A dag object
#' @param result Optional mi_test_result object
#'
#' @return A list of ggplot objects or plots them directly
#'
#' @export
plot_diagnostics <- function(data, dag, result = NULL) {
  data <- validate_data(data, dag)
  n_vars <- ncol(data)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R diagnostics
    par(mfrow = c(2, 2))

    # Marginal distributions
    for (j in 1:min(4, n_vars)) {
      hist(data[, j], main = paste("Marginal:", dag$nodes[j]),
           xlab = dag$nodes[j], col = "lightblue")
    }

    par(mfrow = c(1, 1))
    return(invisible(NULL))
  }

  plots <- list()

  # 1. Marginal distributions
  df_long <- data.frame(
    value = as.vector(data),
    variable = rep(dag$nodes, each = nrow(data))
  )

  plots$marginals <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::labs(title = "Marginal Distributions", x = "Value", y = "Count") +
    ggplot2::theme_minimal()

  # 2. Pairwise scatter plots (first few pairs)
  if (n_vars >= 2) {
    pairs_to_plot <- min(6, n_vars * (n_vars - 1) / 2)
    pair_list <- list()
    count <- 0

    for (i in 1:(n_vars - 1)) {
      for (j in (i + 1):n_vars) {
        count <- count + 1
        if (count > pairs_to_plot) break

        pair_list[[count]] <- data.frame(
          x = data[, i],
          y = data[, j],
          xvar = dag$nodes[i],
          yvar = dag$nodes[j]
        )
      }
      if (count > pairs_to_plot) break
    }

    df_pairs <- do.call(rbind, pair_list)
    df_pairs$pair <- paste(df_pairs$xvar, "vs", df_pairs$yvar)

    plots$pairs <- ggplot2::ggplot(df_pairs, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(alpha = 0.5, color = "steelblue") +
      ggplot2::facet_wrap(~ pair, scales = "free") +
      ggplot2::labs(title = "Pairwise Relationships", x = "", y = "") +
      ggplot2::theme_minimal()
  }

  # 3. Q-Q plots for normality
  plots$qq <- ggplot2::ggplot(df_long, ggplot2::aes(sample = .data$value)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = "red") +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::labs(title = "Q-Q Plots", x = "Theoretical", y = "Sample") +
    ggplot2::theme_minimal()

  # Print or return
  if (interactive()) {
    for (p in plots) print(p)
  }

  invisible(plots)
}
