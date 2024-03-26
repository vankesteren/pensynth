test_that("Data simulation works", {
  dat <- simulate_data(
    N_donor = 20,
    N_covar = 5,
    N_pre = 4,
    N_post = 6
  )
  expect_equal(length(dat$w), 20)
  expect_equal(dim(dat$X0), c(5, 20))
  expect_equal(dim(dat$X1), c(5, 1))
  expect_equal(dim(dat$Y0), c(6, 20))
  expect_equal(dim(dat$Y1), c(6, 1))
  expect_equal(dim(dat$Z0), c(4, 20))
  expect_equal(dim(dat$Z1), c(4, 1))
})
