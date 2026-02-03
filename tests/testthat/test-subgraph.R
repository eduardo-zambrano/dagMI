# Tests for subgraph functions

test_that("extract_subgraph() extracts correct subgraph", {
  A <- matrix(c(0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
                0, 0, 0, 0), 4, 4, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))

  sub_g <- extract_subgraph(g, c("X1", "X2", "X3"))

  expect_equal(sub_g$n, 3)
  expect_equal(sub_g$nodes, c("X1", "X2", "X3"))
  expect_equal(sum(sub_g$adjacency), 2)  # X1->X2, X2->X3
})

test_that("extract_ancestral_subgraph() includes all ancestors", {
  A <- matrix(c(0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
                0, 0, 0, 0), 4, 4, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))

  anc_g <- extract_ancestral_subgraph(g, "X4")

  # Should include all nodes since X4 has X1, X2, X3 as ancestors
  expect_equal(anc_g$n, 4)
})

test_that("extract_subgraph() validates node names", {
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  expect_error(extract_subgraph(g, c("X1", "X4")), "not in DAG")
})

test_that("subgraph_test() tests multiple subgraphs", {
  skip_on_cran()

  set.seed(123)
  A <- matrix(c(0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
                0, 0, 0, 0), 4, 4, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))

  data <- matrix(rnorm(400), 100, 4)
  colnames(data) <- c("X1", "X2", "X3", "X4")

  result <- subgraph_test(data, g, k = 3, B = 20, verbose = FALSE)

  expect_s3_class(result, "subgraph_test_result")
  expect_equal(length(result$p_values), result$n_subgraphs)
  expect_equal(length(result$adjusted_p), result$n_subgraphs)
})

test_that("subgraph_test() applies corrections correctly", {
  skip_on_cran()

  set.seed(123)
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  data <- matrix(rnorm(300), 100, 3)
  colnames(data) <- c("X1", "X2", "X3")

  # Just test that different corrections run without error
  result_bonf <- subgraph_test(data, g, subgraph_nodes = list(c("X1", "X2", "X3")),
                                 correction = "bonferroni", B = 10, verbose = FALSE)
  result_holm <- subgraph_test(data, g, subgraph_nodes = list(c("X1", "X2", "X3")),
                                 correction = "holm", B = 10, verbose = FALSE)
  result_bh <- subgraph_test(data, g, subgraph_nodes = list(c("X1", "X2", "X3")),
                               correction = "bh", B = 10, verbose = FALSE)

  expect_equal(result_bonf$correction, "bonferroni")
  expect_equal(result_holm$correction, "holm")
  expect_equal(result_bh$correction, "bh")
})

test_that("problematic_subgraphs() returns rejected subgraphs", {
  # Create mock result
  mock_result <- structure(
    list(
      subgraph_nodes = list(c("X1", "X2"), c("X2", "X3"), c("X1", "X3")),
      rejected = c(TRUE, FALSE, TRUE)
    ),
    class = "subgraph_test_result"
  )

  prob <- problematic_subgraphs(mock_result)

  expect_equal(length(prob), 2)
  expect_equal(prob[[1]], c("X1", "X2"))
  expect_equal(prob[[2]], c("X1", "X3"))
})

test_that("problematic_nodes() identifies frequently rejected nodes", {
  # Create mock result
  mock_result <- structure(
    list(
      subgraph_nodes = list(c("X1", "X2"), c("X1", "X3"), c("X2", "X3")),
      rejected = c(TRUE, TRUE, FALSE)
    ),
    class = "subgraph_test_result"
  )

  prob_nodes <- problematic_nodes(mock_result)

  expect_s3_class(prob_nodes, "data.frame")
  expect_true("X1" %in% prob_nodes$node)
  # X1 should have count 2 (in both rejected subgraphs)
  expect_equal(prob_nodes$count[prob_nodes$node == "X1"], 2)
})
