# Tests for DAG class

test_that("dag() creates valid DAG from adjacency matrix", {
  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                0, 0, 0), 3, 3, byrow = TRUE)

  g <- dag(A, nodes = c("X1", "X2", "X3"))

  expect_s3_class(g, "dag")
  expect_equal(g$n, 3)
  expect_equal(g$nodes, c("X1", "X2", "X3"))
  expect_true(is_dag(g))
})

test_that("dag() creates valid DAG from edge list", {
  edges <- data.frame(
    from = c("X1", "X2"),
    to = c("X2", "X3")
  )

  g <- dag(edges = edges)

  expect_s3_class(g, "dag")
  expect_equal(g$n, 3)
})

test_that("dag() detects cycles", {
  # Cycle: X1 -> X2 -> X3 -> X1
  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0), 3, 3, byrow = TRUE)

  expect_error(dag(A), "cycles")
})

test_that("dag() detects self-loops", {
  A <- matrix(c(1, 1, 0,
                0, 0, 1,
                0, 0, 0), 3, 3, byrow = TRUE)

  expect_error(dag(A), "self-loops")
})

test_that("parents() returns correct parents", {
  A <- matrix(c(0, 1, 1,
                0, 0, 1,
                0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  expect_equal(parents(g, "X1"), character(0))
  expect_equal(parents(g, "X2"), "X1")
  expect_equal(sort(parents(g, "X3")), c("X1", "X2"))
})

test_that("children() returns correct children", {
  A <- matrix(c(0, 1, 1,
                0, 0, 1,
                0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  expect_equal(sort(children(g, "X1")), c("X2", "X3"))
  expect_equal(children(g, "X2"), "X3")
  expect_equal(children(g, "X3"), character(0))
})

test_that("ancestors() returns correct ancestors", {
  A <- matrix(c(0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
                0, 0, 0, 0), 4, 4, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))

  expect_equal(ancestors(g, "X1"), character(0))
  expect_equal(ancestors(g, "X2"), "X1")
  expect_equal(sort(ancestors(g, "X4")), c("X1", "X2", "X3"))
})

test_that("descendants() returns correct descendants", {
  A <- matrix(c(0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
                0, 0, 0, 0), 4, 4, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3", "X4"))

  expect_equal(sort(descendants(g, "X1")), c("X2", "X3", "X4"))
  expect_equal(sort(descendants(g, "X2")), c("X3", "X4"))
  expect_equal(descendants(g, "X4"), character(0))
})

test_that("topological_orders() finds all orderings", {
  # Fork: X1 -> X2, X1 -> X3 (two orderings)
  A <- matrix(c(0, 1, 1,
                0, 0, 0,
                0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  orders <- topological_orders(g)

  expect_length(orders, 2)

  # Both should start with X1
  expect_true(all(sapply(orders, `[`, 1) == "X1"))
})

test_that("topological_orders() finds unique ordering for chain", {
  # Chain: X1 -> X2 -> X3
  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  orders <- topological_orders(g)

  expect_length(orders, 1)
  expect_equal(orders[[1]], c("X1", "X2", "X3"))
})

test_that("is_valid_ordering() correctly validates orderings", {
  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  expect_true(is_valid_ordering(g, c("X1", "X2", "X3")))
  expect_false(is_valid_ordering(g, c("X3", "X2", "X1")))
  expect_false(is_valid_ordering(g, c("X1", "X3", "X2")))
})

test_that("print.dag() works", {
  A <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)
  g <- dag(A, nodes = c("X1", "X2", "X3"))

  expect_output(print(g), "DAG with 3 nodes")
  expect_output(print(g), "X1 -> X2")
})
