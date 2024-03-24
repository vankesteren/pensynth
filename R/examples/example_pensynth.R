# generate some data
X0 <- matrix(
  c(1, 1.3,
    0.5, 1.8,
    1.1, 2.4,
    1.8, 1.8,
    1.3, 1.8), 2)
X1 <- matrix(c(0.8, 1.65), 2)

# run classic synthetic control (no penalization)
res <- pensynth(X1, X0)

# plot donor units in covariate space
plot(t(X0), asp = 1, xlab = "X1", ylab = "X2",
     main = "Covariate space plot")
# add the treated unit
points(t(X1), pch = 2)
# add the synthetic control
points(t(X0%*%res$w), pch = 3)

# run synthetic control with penalty
res <- pensynth(X1, X0, lambda = 0.5)
# the resulting synthetic control is
# biased towards its closest neighbours
points(t(X0 %*% res$w), pch = 4)
