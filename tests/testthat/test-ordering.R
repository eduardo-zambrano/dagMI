# Tests for ordering selection functions

test_that("get_topological_orders() returns valid orderings", {
  A <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  orders <- get_topological_orders(g)

  # Check all are valid
  for (ord in orders) {
    expect_true(is_valid_ordering(g, ord))
  }
})

test_that("select_ordering() with 'first' returns first ordering", {
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result <- select_ordering(data, g, method = "first")

  expect_equal(result$method, "first")
  expect_true(is_valid_ordering(g, result$ordering))
})

test_that("select_ordering() with 'optimal' finds max Q_n", {
  set.seed(123)
  A <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result <- select_ordering(data, g, method = "optimal")

  expect_equal(result$method, "optimal")
  expect_equal(unname(result$qn), unname(max(result$all_qn)))
})

test_that("bonferroni_orderings() adjusts p-values correctly", {
  p_values <- c(0.01, 0.03, 0.05)

  result <- bonferroni_orderings(p_values, alpha = 0.05)

  expect_equal(result$adjusted_p, p_values * 3)
  expect_equal(result$min_p, 0.01)
  expect_equal(result$adjusted_min_p, 0.03)
})

test_that("fisher_combine() combines p-values", {
  p_values <- c(0.01, 0.02, 0.03)

  result <- fisher_combine(p_values)

  expect_true(result$combined_p >= 0 && result$combined_p <= 1)
  expect_equal(result$df, 6)  # 2 * 3
  expect_true(result$chi_sq > 0)
})

test_that("stouffer_combine() combines p-values", {
  p_values <- c(0.01, 0.02, 0.03)

  result <- stouffer_combine(p_values)

  expect_true(result$combined_p >= 0 && result$combined_p <= 1)
  expect_true(!is.na(result$z_score))
})

test_that("count_orderings() counts correctly for small DAGs", {
  # Chain: 1 ordering
  A1 <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g1 <- dag(A1, nodes = c("X1", "X2", "X3"))

  expect_equal(count_orderings(g1), 1)

  # Fork: 2 orderings
  A2 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
  g2 <- dag(A2, nodes = c("X1", "X2", "X3"))

  expect_equal(count_orderings(g2), 2)

  # Empty graph with 3 nodes: 6 orderings (3!)
  A3 <- matrix(0, 3, 3)
  g3 <- dag(A3, nodes = c("X1", "X2", "X3"))

  expect_equal(count_orderings(g3), 6)
})
