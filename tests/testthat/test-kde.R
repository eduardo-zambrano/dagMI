# Tests for KDE functions

test_that("kde_marginal() produces valid density", {
  set.seed(123)
  x <- rnorm(100)

  kde <- kde_marginal(x)

  expect_s3_class(kde, "kde_marginal")
  expect_true(all(kde$y >= 0))
  expect_equal(kde$n, 100)
  expect_true(kde$bandwidth > 0)

  # Integral should be approximately 1
  integral <- sum(kde$y) * diff(kde$x[1:2])
  expect_true(abs(integral - 1) < 0.1)
})

test_that("kde_eval() evaluates at specified points", {
  set.seed(123)
  x <- rnorm(100)
  kde <- kde_marginal(x)

  # Evaluate at specific points
  x_new <- c(-1, 0, 1)
  dens <- kde_eval(kde, x_new)

  expect_length(dens, 3)
  expect_true(all(dens > 0))
})

test_that("kde_bivariate() produces valid 2D density", {
  set.seed(123)
  x <- rnorm(100)
  y <- x + rnorm(100, sd = 0.5)

  kde <- kde_bivariate(x, y, n_grid = 32)

  expect_s3_class(kde, "kde_bivariate")
  expect_true(all(kde$z >= 0))
  expect_equal(dim(kde$z), c(32, 32))
})

test_that("kde_conditional() creates conditional KDE object", {
  set.seed(123)
  x <- rnorm(100)
  y <- 2 * x + rnorm(100, sd = 0.5)

  kde <- kde_conditional(y, x)

  expect_s3_class(kde, "kde_conditional")
  expect_equal(kde$n, 100)
})

test_that("select_bandwidth() returns positive bandwidth", {
  set.seed(123)
  x <- rnorm(100)

  bw_silverman <- select_bandwidth(x, "silverman")
  bw_scott <- select_bandwidth(x, "scott")
  bw_cv <- select_bandwidth(x, "cv", n_grid = 10)

  expect_true(bw_silverman > 0)
  expect_true(bw_scott > 0)
  expect_true(bw_cv > 0)
})

test_that("KDE handles edge cases", {
  # Very small sample - still works
  set.seed(123)
  x <- rnorm(5)
  kde <- kde_marginal(x)
  expect_true(all(kde$y >= 0))

  # Data with outliers
  x <- c(rnorm(100), 10, -10)
  kde <- kde_marginal(x)
  expect_true(all(kde$y >= 0))
})
