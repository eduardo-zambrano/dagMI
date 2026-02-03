#' Utility Functions for dagMI Package
#'
#' @name utils
#' @keywords internal
NULL

#' Null coalescing operator
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Validate Data Matrix
#'
#' @param data A numeric matrix or data frame
#' @param dag Optional dag object to check column matching
#' @return The data as a matrix with column names
#' @keywords internal
validate_data <- function(data, dag = NULL) {
  data <- as.matrix(data)

  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }

  if (any(is.na(data))) {
    warning("Data contains NA values, which may affect results")
  }

  if (nrow(data) < 10) {
    warning("Very small sample size (n < 10)")
  }

  if (!is.null(dag)) {
    if (ncol(data) != dag$n) {
      stop("Number of columns in data (", ncol(data),
           ") does not match number of nodes in DAG (", dag$n, ")")
    }

    if (!is.null(colnames(data))) {
      if (!all(colnames(data) == dag$nodes)) {
        # Try to reorder
        if (all(dag$nodes %in% colnames(data))) {
          data <- data[, dag$nodes, drop = FALSE]
        } else {
          warning("Column names do not match DAG nodes; assuming same order")
        }
      }
    } else {
      colnames(data) <- dag$nodes
    }
  }

  data
}

#' Silverman's Rule of Thumb for Bandwidth
#'
#' @param x Numeric vector
#' @return Bandwidth value
#' @keywords internal
silverman_bandwidth <- function(x) {
  n <- length(x)
  sigma <- sd(x)
  iqr <- IQR(x)

  # Silverman's rule
  h <- 0.9 * min(sigma, iqr / 1.34) * n^(-1/5)

  # Ensure positive
  max(h, 1e-6)
}

#' Scott's Rule for Bandwidth
#'
#' @param x Numeric vector
#' @return Bandwidth value
#' @keywords internal
scott_bandwidth <- function(x) {
  n <- length(x)
  sigma <- sd(x)
  h <- sigma * n^(-1/5)
  max(h, 1e-6)
}

#' Gauss-Legendre Quadrature Nodes and Weights
#'
#' @param n Number of quadrature points
#' @return List with nodes and weights
#' @keywords internal
gauss_legendre <- function(n) {
  # Use Golub-Welsch algorithm or lookup for small n
  if (n <= 20) {
    return(.gauss_legendre_lookup(n))
  }

  # For larger n, use iterative algorithm
  .gauss_legendre_compute(n)
}

#' Lookup Table for Small Gauss-Legendre Rules
#' @keywords internal
.gauss_legendre_lookup <- function(n) {
  # Precomputed values for common orders
  rules <- list(
    `1` = list(
      nodes = 0,
      weights = 2
    ),
    `2` = list(
      nodes = c(-0.5773502691896257, 0.5773502691896257),
      weights = c(1, 1)
    ),
    `3` = list(
      nodes = c(-0.7745966692414834, 0, 0.7745966692414834),
      weights = c(0.5555555555555556, 0.8888888888888888, 0.5555555555555556)
    ),
    `4` = list(
      nodes = c(-0.8611363115940526, -0.3399810435848563,
                0.3399810435848563, 0.8611363115940526),
      weights = c(0.3478548451374538, 0.6521451548625461,
                  0.6521451548625461, 0.3478548451374538)
    ),
    `5` = list(
      nodes = c(-0.9061798459386640, -0.5384693101056831, 0,
                0.5384693101056831, 0.9061798459386640),
      weights = c(0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                  0.4786286704993665, 0.2369268850561891)
    ),
    `10` = list(
      nodes = c(-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
                0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
                0.9739065285171717),
      weights = c(0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
                  0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
                  0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
                  0.0666713443086881)
    )
  )

  key <- as.character(n)
  if (key %in% names(rules)) {
    return(rules[[key]])
  }

  # Fall back to computation
  .gauss_legendre_compute(n)
}

#' Compute Gauss-Legendre Nodes and Weights
#' @keywords internal
.gauss_legendre_compute <- function(n) {
  # Newton's method to find roots of Legendre polynomials
  nodes <- numeric(n)
  weights <- numeric(n)

  m <- (n + 1) %/% 2

  for (i in 1:m) {
    # Initial guess
    z <- cos(pi * (i - 0.25) / (n + 0.5))

    # Newton iteration
    for (iter in 1:100) {
      p1 <- 1
      p2 <- 0

      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
      }

      pp <- n * (z * p1 - p2) / (z^2 - 1)
      z_new <- z - p1 / pp

      if (abs(z_new - z) < 1e-15) break
      z <- z_new
    }

    nodes[i] <- -z
    nodes[n + 1 - i] <- z
    weights[i] <- 2 / ((1 - z^2) * pp^2)
    weights[n + 1 - i] <- weights[i]
  }

  list(nodes = nodes, weights = weights)
}

#' Transform Quadrature to Arbitrary Interval
#'
#' @param gl Gauss-Legendre rule (list with nodes and weights)
#' @param a Lower bound
#' @param b Upper bound
#' @return Transformed rule
#' @keywords internal
transform_quadrature <- function(gl, a, b) {
  mid <- (a + b) / 2
  half_width <- (b - a) / 2

  list(
    nodes = mid + half_width * gl$nodes,
    weights = half_width * gl$weights
  )
}

#' Create Tensor Product Quadrature Grid
#'
#' @param n_points Number of points per dimension
#' @param dims Number of dimensions
#' @param bounds Matrix with 2 columns (lower, upper) and dims rows
#' @return List with grid matrix and weights vector
#' @keywords internal
tensor_quadrature <- function(n_points, dims, bounds) {
  gl <- gauss_legendre(n_points)

  # Transform to each dimension's bounds
  rules <- lapply(1:dims, function(d) {
    transform_quadrature(gl, bounds[d, 1], bounds[d, 2])
  })

  # Create tensor product
  grids <- lapply(rules, function(r) r$nodes)
  wts <- lapply(rules, function(r) r$weights)

  grid <- as.matrix(expand.grid(grids))
  weight_grid <- expand.grid(wts)
  weights <- apply(weight_grid, 1, prod)

  list(grid = grid, weights = weights)
}

#' Safe Log for Small Values
#'
#' @param x Numeric vector
#' @param min_val Minimum value before log
#' @return log(max(x, min_val))
#' @keywords internal
safe_log <- function(x, min_val = 1e-300) {
  log(pmax(x, min_val))
}

#' Safe Division
#'
#' @param num Numerator
#' @param denom Denominator
#' @param default Value when denominator is near zero
#' @return num/denom or default
#' @keywords internal
safe_divide <- function(num, denom, default = 0) {
  result <- num / denom
  result[abs(denom) < 1e-15] <- default
  result
}

#' Progress Bar Wrapper
#'
#' @param total Total iterations
#' @param verbose Whether to show progress
#' @return Progress bar object or NULL
#' @keywords internal
make_progress <- function(total, verbose = TRUE) {
  if (!verbose) return(NULL)

  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps = total)
  } else {
    # Simple text progress
    list(
      .env = new.env(),
      update = function() {
        if (!exists("count", envir = .env)) .env$count <- 0
        .env$count <- .env$count + 1
        if (.env$count %% ceiling(total / 20) == 0) {
          cat(sprintf("\r  Progress: %d%%", round(100 * .env$count / total)))
        }
      },
      finish = function() cat("\n")
    )
  }
}

#' Set Up Parallel Backend
#'
#' @param parallel Logical or number of cores
#' @return Previous plan (for restoration)
#' @keywords internal
setup_parallel <- function(parallel = TRUE) {
  if (!requireNamespace("future", quietly = TRUE)) {
    warning("Package 'future' not available; running sequentially")
    return(NULL)
  }

  old_plan <- future::plan()

  if (isTRUE(parallel)) {
    n_cores <- parallel::detectCores() - 1
    n_cores <- max(1, n_cores)
  } else if (is.numeric(parallel) && parallel > 1) {
    n_cores <- parallel
  } else {
    return(NULL)
  }

  future::plan(future::multisession, workers = n_cores)
  old_plan
}

#' Restore Parallel Backend
#'
#' @param old_plan Previous plan from setup_parallel
#' @keywords internal
restore_parallel <- function(old_plan) {
  if (!is.null(old_plan) && requireNamespace("future", quietly = TRUE)) {
    future::plan(old_plan)
  }
}

#' Check if Package is Available
#'
#' @param pkg Package name
#' @param purpose What the package is needed for
#' @return Logical
#' @keywords internal
check_package <- function(pkg, purpose = NULL) {
  available <- requireNamespace(pkg, quietly = TRUE)
  if (!available && !is.null(purpose)) {
    message("Install '", pkg, "' for ", purpose)
  }
  available
}
