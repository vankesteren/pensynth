test_that("convhull detection works", {
  set.seed(45)
  N  <- 1000
  P  <- 9
  X0 <- matrix(runif(N*P), P)
  X1 <- matrix(rep(0.5, P))
  expect_true(in_convex_hull(X1, X0))
  X1 <- matrix(rep(.99, P))
  expect_false(in_convex_hull(X1, X0))
})
