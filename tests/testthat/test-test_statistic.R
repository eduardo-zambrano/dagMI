# Tests for test statistic computation

test_that("h_squared() computes x^2", {
  x <- c(-2, -1, 0, 1, 2)
  expect_equal(h_squared(x), c(4, 1, 0, 1, 4))
})

test_that("h_poly() computes x^(2k)", {
  x <- c(-2, -1, 0, 1, 2)
  expect_equal(h_poly(x, k = 1), x^2)
  expect_equal(h_poly(x, k = 2), x^4)
})

test_that("h_exp() computes exp(t*x)", {
  x <- c(-1, 0, 1)
  expect_equal(h_exp(x, t = 1), exp(x))
  expect_equal(h_exp(x, t = 0.5), exp(0.5 * x))
})

test_that("h_indicator() computes I(x > c)", {
  x <- c(-1, 0, 0.5, 1, 2)
  expect_equal(h_indicator(x, c = 0), c(0, 0, 1, 1, 1))
  expect_equal(h_indicator(x, c = 1), c(0, 0, 0, 0, 1))
})

test_that("compute_test_stat() returns valid result", {
  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  X1 <- rnorm(100)
  X2 <- 0.5 * X1 + rnorm(100, sd = 0.5)
  X3 <- 0.5 * X2 + rnorm(100, sd = 0.5)
  data <- cbind(X1, X2, X3)

  result <- compute_test_stat(data, g)

  expect_s3_class(result, "test_stat_result")
  expect_true(is.numeric(result$statistic))
  expect_true(is.numeric(result$lhs))
  expect_true(is.numeric(result$rhs))
  expect_equal(result$h_name, "squared")
})

test_that("compute_test_stat() works with different test functions", {
  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result_sq <- compute_test_stat(data, g, h = "squared")
  result_poly <- compute_test_stat(data, g, h = "poly", h_params = list(k = 2))
  result_exp <- compute_test_stat(data, g, h = "exp", h_params = list(t = 0.5))

  expect_equal(result_sq$h_name, "squared")
  expect_equal(result_poly$h_name, "poly")
  expect_equal(result_exp$h_name, "exp")
})

test_that("compute_test_stat() works with custom function", {
  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  custom_h <- function(x) abs(x)^1.5

  result <- compute_test_stat(data, g, h = custom_h)

  expect_equal(result$h_name, "custom")
})

test_that("check_admissibility() correctly identifies admissible functions", {
  set.seed(123)
  data <- matrix(rnorm(300), 100, 3)

  # h_squared is admissible
  adm <- check_admissibility(h_squared, data)
  expect_true(adm$admissible)
  expect_true(adm$non_negative)
  expect_true(adm$finite_moments)

  # identity function (can be negative) is not admissible
  adm_neg <- check_admissibility(identity, data)
  expect_false(adm_neg$admissible)
  expect_false(adm_neg$non_negative)
})
