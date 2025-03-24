set.seed(45)
N_covar <- 7
N_donor <- 50
N_target <- 12

w  <- runif(N_donor)
w[5:N_donor] <- 0
w  <- w / sum(w)
v  <- rep(1, N_covar)
X0 <- matrix(rnorm(N_covar*N_donor), N_covar)
X1 <- X0%*%w
Z0 <- matrix(rnorm(N_target*N_donor), N_target)
Z1 <- Z0%*%w

test_that("CV pensynth works", {
  res <- cv_pensynth(X1, X0, Z1, Z0, v, verbose = FALSE)
  expect_lt(crossprod(res$w_opt - w), 5e-3)
  expect_gt(max(res$w_path[,100]), 0.999)
  expect_lt(max(res$w_path[,1]), 0.5)
})

