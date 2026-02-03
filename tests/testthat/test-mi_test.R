# Tests for main mi_test() function

test_that("mi_test() returns valid result with bootstrap", {
  skip_on_cran()  # Skip slow test on CRAN

  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  X1 <- rnorm(100)
  X2 <- 0.5 * X1 + rnorm(100, sd = 0.5)
  X3 <- 0.5 * X2 + rnorm(100, sd = 0.5)
  data <- cbind(X1, X2, X3)

  result <- mi_test(data, g, B = 50, verbose = FALSE)

  expect_s3_class(result, "mi_test_result")
  expect_true(!is.na(result$p_value))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$decision %in% c("reject", "fail_to_reject"))
})

test_that("mi_test() works without bootstrap (B = 0)", {
  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result <- mi_test(data, g, B = 0, verbose = FALSE)

  expect_true(is.na(result$p_value))
  expect_equal(result$decision, "not_computed")
})

test_that("mi_test() respects ordering parameter", {
  set.seed(123)
  A <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  result1 <- mi_test(data, g, ordering = c("X1", "X2", "X3"),
                      B = 0, verbose = FALSE)
  result2 <- mi_test(data, g, ordering = c("X1", "X3", "X2"),
                      B = 0, verbose = FALSE)

  expect_equal(result1$ordering, c("X1", "X2", "X3"))
  expect_equal(result2$ordering, c("X1", "X3", "X2"))
})

test_that("compare_dags() compares two DAGs", {
  skip_on_cran()

  set.seed(123)
  # Chain DAG
  A1 <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g1 <- dag(A1, nodes = c("X1", "X2", "X3"))

  # Fork DAG
  A2 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)
  g2 <- dag(A2, nodes = c("X1", "X2", "X3"))

  # Generate data from chain
  X1 <- rnorm(100)
  X2 <- 0.5 * X1 + rnorm(100, sd = 0.5)
  X3 <- 0.5 * X2 + rnorm(100, sd = 0.5)
  data <- cbind(X1, X2, X3)

  result <- compare_dags(data, g1, g2, B = 50, verbose = FALSE)

  expect_s3_class(result, "dag_comparison")
  expect_true(!is.null(result$comparison$preferred))
})

test_that("mi_test() handles data validation", {
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  # Wrong number of columns
  data <- matrix(rnorm(200), 100, 2)
  expect_error(mi_test(data, g, B = 0, verbose = FALSE), "columns")

  # Non-numeric data
  data_char <- data.frame(a = letters[1:10], b = letters[11:20], c = letters[1:10])
  expect_error(mi_test(data_char, g, B = 0, verbose = FALSE), "numeric")
})
