#' Kernel Density Estimation Functions
#'
#' @description
#' Functions for univariate, bivariate, and conditional kernel density
#' estimation with Gaussian kernels.
#'
#' @name kde
NULL

#' Univariate Kernel Density Estimation
#'
#' @description
#' Estimates a univariate density using Gaussian kernel.
#'
#' @param x Numeric vector of data points
#' @param bandwidth Bandwidth parameter. If NULL, uses Silverman's rule.
#' @param n_grid Number of grid points for density evaluation (default 512)
#' @param from Lower bound of evaluation range (default: min(x) - 3*bw)
#' @param to Upper bound of evaluation range (default: max(x) + 3*bw)
#'
#' @return An object of class "kde_marginal" with components:
#'   \item{x}{Grid points}
#'   \item{y}{Density values}
#'   \item{data}{Original data}
#'   \item{bandwidth}{Bandwidth used}
#'   \item{n}{Sample size}
#'
#' @details
#' Uses the Gaussian kernel K(u) = (2*pi)^(-1/2) * exp(-u^2/2).
#' The default bandwidth is Silverman's rule of thumb.
#'
#' @examples
#' x <- rnorm(100)
#' kde <- kde_marginal(x)
#' plot(kde$x, kde$y, type = "l")
#'
#' @export
kde_marginal <- function(x, bandwidth = NULL, n_grid = 512,
                          from = NULL, to = NULL) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 2) {
    stop("Need at least 2 data points")
  }

  # Bandwidth selection
  if (is.null(bandwidth)) {
    bandwidth <- silverman_bandwidth(x)
  }

  # Set evaluation range
  if (is.null(from)) from <- min(x) - 3 * bandwidth
  if (is.null(to)) to <- max(x) + 3 * bandwidth

  # Create evaluation grid
  x_grid <- seq(from, to, length.out = n_grid)

  # Evaluate density using C++
  y <- kde_eval_cpp(x_grid, x, bandwidth)

  structure(
    list(
      x = x_grid,
      y = y,
      data = x,
      bandwidth = bandwidth,
      n = n
    ),
    class = "kde_marginal"
  )
}

#' Evaluate KDE at Specific Points
#'
#' @description
#' Evaluates a fitted KDE at arbitrary points.
#'
#' @param kde A kde_marginal object
#' @param x_new Points at which to evaluate density
#' @param truncate Minimum density value (default 1e-10)
#'
#' @return Numeric vector of density values
#'
#' @examples
#' kde <- kde_marginal(rnorm(100))
#' kde_eval(kde, c(-1, 0, 1))
#'
#' @export
kde_eval <- function(kde, x_new, truncate = 1e-10) {
  if (!inherits(kde, "kde_marginal")) {
    stop("'kde' must be a kde_marginal object")
  }

  density <- kde_eval_truncated_cpp(x_new, kde$data, kde$bandwidth, truncate)
  as.numeric(density)
}

#' Bivariate Kernel Density Estimation
#'
#' @description
#' Estimates a bivariate density using product Gaussian kernel.
#'
#' @param x Numeric vector for first variable
#' @param y Numeric vector for second variable
#' @param bandwidth Length-2 vector of bandwidths. If NULL, uses Silverman's rule
#'   for each variable.
#' @param n_grid Number of grid points per dimension (default 64)
#'
#' @return An object of class "kde_bivariate" with components:
#'   \item{x}{Grid for first variable}
#'   \item{y}{Grid for second variable}
#'   \item{z}{Matrix of density values}
#'   \item{data}{Original data matrix}
#'   \item{bandwidth}{Bandwidths used}
#'   \item{n}{Sample size}
#'
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100, sd = 0.5)
#' kde <- kde_bivariate(x, y)
#'
#' @export
kde_bivariate <- function(x, y, bandwidth = NULL, n_grid = 64) {
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have same length")
  }

  # Remove NAs
  complete <- !is.na(x) & !is.na(y)
  x <- x[complete]
  y <- y[complete]
  n <- length(x)

  if (n < 2) {
    stop("Need at least 2 complete observations")
  }

  # Bandwidth selection
  if (is.null(bandwidth)) {
    bandwidth <- c(silverman_bandwidth(x), silverman_bandwidth(y))
  } else if (length(bandwidth) == 1) {
    bandwidth <- rep(bandwidth, 2)
  }

  # Create grid
  x_range <- range(x) + c(-3, 3) * bandwidth[1]
  y_range <- range(y) + c(-3, 3) * bandwidth[2]

  x_grid <- seq(x_range[1], x_range[2], length.out = n_grid)
  y_grid <- seq(y_range[1], y_range[2], length.out = n_grid)

  # Create evaluation points
  eval_grid <- expand.grid(x = x_grid, y = y_grid)
  eval_matrix <- as.matrix(eval_grid)
  data_matrix <- cbind(x, y)

  # Evaluate density
  density <- kde_bivariate_eval_cpp(eval_matrix, data_matrix, bandwidth)

  z <- matrix(density, nrow = n_grid, ncol = n_grid)

  structure(
    list(
      x = x_grid,
      y = y_grid,
      z = z,
      data = data_matrix,
      bandwidth = bandwidth,
      n = n
    ),
    class = "kde_bivariate"
  )
}

#' Conditional Density Estimation
#'
#' @description
#' Estimates conditional density p(y|x) using kernel methods.
#'
#' @param y Response variable (what we want the density of)
#' @param x Conditioning variable
#' @param bandwidth_y Bandwidth for y. If NULL, uses Silverman's rule.
#' @param bandwidth_x Bandwidth for x. If NULL, uses Silverman's rule.
#' @param n_grid Number of grid points for y (default 512)
#'
#' @return An object of class "kde_conditional" with components:
#'   \item{y_grid}{Grid for y}
#'   \item{data_y}{Original y data}
#'   \item{data_x}{Original x data}
#'   \item{bandwidth_y}{Bandwidth for y}
#'   \item{bandwidth_x}{Bandwidth for x}
#'   \item{n}{Sample size}
#'
#' @details
#' Uses the Nadaraya-Watson estimator approach where
#' p(y|x) = f(x,y) / f(x) is estimated via kernel smoothing.
#'
#' @examples
#' x <- rnorm(100)
#' y <- 2*x + rnorm(100, sd = 0.5)
#' kde <- kde_conditional(y, x)
#'
#' @export
kde_conditional <- function(y, x, bandwidth_y = NULL, bandwidth_x = NULL,
                             n_grid = 512) {
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have same length")
  }

  # Remove NAs
  complete <- !is.na(x) & !is.na(y)
  x <- x[complete]
  y <- y[complete]
  n <- length(x)

  if (n < 2) {
    stop("Need at least 2 complete observations")
  }

  # Bandwidth selection
  if (is.null(bandwidth_y)) bandwidth_y <- silverman_bandwidth(y)
  if (is.null(bandwidth_x)) bandwidth_x <- silverman_bandwidth(x)

  # Create y grid
  y_range <- range(y) + c(-3, 3) * bandwidth_y
  y_grid <- seq(y_range[1], y_range[2], length.out = n_grid)

  structure(
    list(
      y_grid = y_grid,
      data_y = y,
      data_x = x,
      bandwidth_y = bandwidth_y,
      bandwidth_x = bandwidth_x,
      n = n
    ),
    class = "kde_conditional"
  )
}

#' Evaluate Conditional Density
#'
#' @description
#' Evaluates conditional density p(y|x) at specified points.
#'
#' @param kde A kde_conditional object
#' @param y_new Points at which to evaluate y density
#' @param x_cond Conditioning values for x
#' @param truncate Minimum density value (default 1e-10)
#'
#' @return Numeric vector of conditional density values
#'
#' @export
kde_cond_eval <- function(kde, y_new, x_cond, truncate = 1e-10) {
  if (!inherits(kde, "kde_conditional")) {
    stop("'kde' must be a kde_conditional object")
  }

  # Ensure same length
  if (length(y_new) != length(x_cond)) {
    if (length(x_cond) == 1) {
      x_cond <- rep(x_cond, length(y_new))
    } else {
      stop("'y_new' and 'x_cond' must have same length")
    }
  }

  density <- kde_conditional_cpp(
    y_new, x_cond,
    kde$data_y, kde$data_x,
    kde$bandwidth_x, kde$bandwidth_y,
    truncate
  )

  as.numeric(density)
}

#' Bandwidth Selection
#'
#' @description
#' Selects bandwidth for kernel density estimation.
#'
#' @param x Numeric vector of data
#' @param method Bandwidth selection method:
#'   \itemize{
#'     \item "silverman": Silverman's rule of thumb (default)
#'     \item "scott": Scott's rule
#'     \item "cv": Leave-one-out cross-validation
#'   }
#' @param n_grid Number of bandwidths to try for CV (default 30)
#'
#' @return Optimal bandwidth value
#'
#' @examples
#' x <- rnorm(100)
#' select_bandwidth(x, "silverman")
#' select_bandwidth(x, "cv")
#'
#' @export
select_bandwidth <- function(x, method = "silverman", n_grid = 30) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]

  method <- match.arg(method, c("silverman", "scott", "cv"))

  if (method == "silverman") {
    return(silverman_bandwidth(x))
  }

  if (method == "scott") {
    return(scott_bandwidth(x))
  }

  # Cross-validation
  h_silverman <- silverman_bandwidth(x)

  # Search range around Silverman
  h_grid <- exp(seq(log(h_silverman / 3), log(h_silverman * 3),
                    length.out = n_grid))

  cv_scores <- sapply(h_grid, function(h) {
    kde_cv_score_cpp(x, h)
  })

  h_grid[which.max(cv_scores)]
}

#' Create KDE Objects for All Variables
#'
#' @description
#' Creates marginal KDE objects for each column of a data matrix.
#'
#' @param data Numeric matrix or data frame
#' @param bandwidth Vector of bandwidths. If NULL, uses Silverman's rule.
#'
#' @return List of kde_marginal objects
#'
#' @keywords internal
create_kde_list <- function(data, bandwidth = NULL) {
  data <- as.matrix(data)
  n_vars <- ncol(data)

  if (is.null(bandwidth)) {
    bandwidth <- apply(data, 2, silverman_bandwidth)
  } else if (length(bandwidth) == 1) {
    bandwidth <- rep(bandwidth, n_vars)
  }

  kdes <- lapply(seq_len(n_vars), function(j) {
    kde_marginal(data[, j], bandwidth = bandwidth[j])
  })

  # Also store as list format for C++
  kde_objects <- lapply(seq_len(n_vars), function(j) {
    list(
      data = data[, j],
      bandwidth = bandwidth[j]
    )
  })

  list(kdes = kdes, kde_objects = kde_objects)
}

#' Sample from Conditional Distribution
#'
#' @description
#' Generates samples from an estimated conditional distribution p(y|x).
#'
#' @param kde A kde_conditional object
#' @param x_cond Conditioning values
#' @param n_samples Number of samples per conditioning value
#'
#' @return Matrix of samples (rows = samples, cols = conditioning values)
#'
#' @details
#' Uses rejection sampling with the marginal density of y as proposal.
#'
#' @export
kde_cond_sample <- function(kde, x_cond, n_samples = 1) {
  if (!inherits(kde, "kde_conditional")) {
    stop("'kde' must be a kde_conditional object")
  }

  n_cond <- length(x_cond)
  samples <- matrix(NA, nrow = n_samples, ncol = n_cond)

  for (j in seq_len(n_cond)) {
    samples[, j] <- .rejection_sample_conditional(
      kde, x_cond[j], n_samples
    )
  }

  samples
}

#' Rejection Sampling from Conditional
#' @keywords internal
.rejection_sample_conditional <- function(kde, x_cond, n_samples) {
  # Use marginal of y as proposal
  y_range <- range(kde$data_y)
  y_range <- y_range + c(-3, 3) * kde$bandwidth_y

  # Find approximate max of conditional density
  y_test <- seq(y_range[1], y_range[2], length.out = 100)
  x_test <- rep(x_cond, 100)
  f_test <- kde_cond_eval(kde, y_test, x_test)
  M <- max(f_test) * 1.5  # Envelope constant

  samples <- numeric(n_samples)
  count <- 0
  max_iter <- n_samples * 100

  iter <- 0
  while (count < n_samples && iter < max_iter) {
    iter <- iter + 1

    # Propose from uniform over range
    y_prop <- runif(1, y_range[1], y_range[2])

    # Evaluate density
    f_prop <- kde_cond_eval(kde, y_prop, x_cond)

    # Accept/reject
    u <- runif(1)
    if (u < f_prop / M) {
      count <- count + 1
      samples[count] <- y_prop
    }
  }

  if (count < n_samples) {
    warning("Rejection sampling did not find enough samples")
    # Fill remaining with jittered data
    samples[(count + 1):n_samples] <- sample(kde$data_y, n_samples - count,
                                              replace = TRUE) +
      rnorm(n_samples - count, sd = kde$bandwidth_y)
  }

  samples
}
