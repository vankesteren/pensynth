set.seed(123)

w <- c(.5, .4, .1, 0, 0, 0)
X0 <- matrix(rnorm(3*6), 3)
X1 <- X0 %*% w


test_that("Basic pensynth works", {
  res <- pensynth(X1, X0, 1, lambda = 0)
  expect_length(res$w, 6)
  expect_lt(sum(abs(res$w - w)), 5e-3)
})

test_that("Pensynth with high lambda works", {
  res <- pensynth(X1, X0, 1, lambda = 100)
  expect_gt(max(res$w), 0.99)
})

