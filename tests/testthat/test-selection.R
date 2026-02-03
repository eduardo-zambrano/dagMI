# Tests for selection.R functions

test_that("select_test_function works correctly", {
  # Create chain DAG
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  chain <- dag(A, nodes = c("X1", "X2", "X3"))

  # Generate data
  set.seed(42)
  n <- 150
  X1 <- rnorm(n)
  X2 <- 0.5 * X1 + rnorm(n, sd = 0.5)
  X3 <- 0.5 * X2 + rnorm(n, sd = 0.5)
  data <- cbind(X1, X2, X3)

  # Test adaptive selection
  result <- select_test_function(data, chain, B = 30, verbose = FALSE, seed = 123)

  expect_s3_class(result, "select_test_result")
  expect_true(result$selected %in% names(result$candidates))
  expect_true(result$n_selection + result$n_inference == n)
  expect_true(!is.na(result$final_result$p_value))
})

test_that("power_mi_test computes power correctly", {
  # Create chain DAG
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  chain <- dag(A, nodes = c("X1", "X2", "X3"))

  # Type I error test (same DAG)
  power_result <- power_mi_test(
    dag_true = chain,
    dag_test = chain,
    n = 80,
    n_sim = 15,
    B = 30,
    verbose = FALSE,
    seed = 456
  )

  expect_s3_class(power_result, "power_result")
  expect_true(power_result$same_dag)
  expect_true(power_result$power >= 0 && power_result$power <= 1)
  expect_length(power_result$p_values, 15)
})

test_that("diagnose_confounding identifies issues", {
  # Create DAGs
  A_chain <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  A_fork <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)

  chain <- dag(A_chain, nodes = c("X1", "X2", "X3"))
  fork <- dag(A_fork, nodes = c("X1", "X2", "X3"))

  # Generate data from chain
  set.seed(789)
  n <- 150
  X1 <- rnorm(n)
  X2 <- 0.5 * X1 + rnorm(n, sd = 0.5)
  X3 <- 0.5 * X2 + rnorm(n, sd = 0.5)
  data <- cbind(X1, X2, X3)

  # Diagnose
  diag <- diagnose_confounding(data, list(chain, fork), B = 30, verbose = FALSE)

  expect_s3_class(diag, "confounding_diagnostic")
  expect_equal(diag$n_dags, 2)
  expect_true(diag$diagnosis %in% c("ALL_REJECTED", "SOME_REJECTED", "NONE_REJECTED"))
  expect_length(diag$recommendations, length(diag$recommendations))
})
