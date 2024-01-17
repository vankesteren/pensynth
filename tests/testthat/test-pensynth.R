X0 <- matrix(
  c(1, 1.3,
    0.5, 1.8,
    1.1, 2.4,
    1.8, 1.8,
    1.3, 1.8), 2)
X1 <- matrix(c(0.8, 1.65), 2)
v <- rep(1, 2)

test_that("Basic pensynth works", {
  res <- pensynth::pensynth(X1, X0, v, lambda = 0)
  expect_length(res$w, 5)
  expect_lt(sum(abs(res$w - c(0.344, 0.504, 0.036, 0.028, 0.087))), 5e-3)
})

test_that("Pensynth with high lambda works", {
  res <- pensynth(X1, X0, v, lambda = 100)
  expect_gt(max(res$w), 0.99)
})

