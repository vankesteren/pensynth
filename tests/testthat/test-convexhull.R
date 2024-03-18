set.seed(45)
N  <- 1000
P  <- 9
X0 <- matrix(runif(N*P), P)

test_that("convhull detection works", {
  X1 <- matrix(rep(0.5, P))
  expect_true(in_convex_hull(X1, X0))
  X1 <- matrix(rep(.99, P))
  expect_false(in_convex_hull(X1, X0))
})

test_that("convhull with multiple treated units works", {
  X1 <- cbind(rep(0.5, P), rep(.99, P))
  expect_equal(in_convex_hull(X1, X0), c(TRUE, FALSE))
})
