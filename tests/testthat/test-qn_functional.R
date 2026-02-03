# Tests for Q_n functional computation

test_that("compute_qn() returns valid qn_result", {
  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  X1 <- rnorm(100)
  X2 <- 0.5 * X1 + rnorm(100, sd = 0.5)
  X3 <- 0.5 * X2 + rnorm(100, sd = 0.5)
  data <- cbind(X1, X2, X3)

  result <- compute_qn(data, g)

  expect_s3_class(result, "qn_result")
  expect_true(result$qn > 0)
  expect_equal(result$n_vars, 3)
  expect_equal(result$ordering, c("X1", "X2", "X3"))
})

test_that("compute_qn() works with different methods", {
  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result_quad <- compute_qn(data, g, method = "quadrature", n_quad = 5)
  result_mc <- compute_qn(data, g, method = "montecarlo", n_mc = 1000)

  expect_equal(result_quad$method, "quadrature")
  expect_equal(result_mc$method, "montecarlo")

  # Monte Carlo should have standard error
  expect_false(is.na(result_mc$se))
})

test_that("compute_qn_all_orderings() evaluates multiple orderings", {
  set.seed(123)
  # Fork: two orderings possible
  A <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result <- compute_qn_all_orderings(data, g)

  expect_equal(length(result$results), 2)
  expect_equal(length(result$qn_values), 2)
  expect_true(!is.null(result$best_ordering))
  expect_equal(unname(result$best_qn), unname(max(result$qn_values)))
})

test_that("qn_bounds() returns valid bounds", {
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  bounds <- qn_bounds(g)

  expect_equal(bounds$lower, 0)
  expect_equal(bounds$upper, 1)
})
