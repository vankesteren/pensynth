# visual exploration of exponential weighting

X0 <- matrix(c(1, 2, 3, 3, 2, 1), ncol = 3)
X1 <- matrix(c(2.5, 1.5))



par(mfrow = c(3, 3))

plot(rbind(t(X0), t(X1)), asp = 1, pch = c("a", "b", "c", "t"), xlab = "x", ylab = "y")
polygon(t(X0), lty = 3)

replicate(8, {
  w0 <- rexp(3)
  w1 <- rexp(1)
  w0 <- w0 / sum(w0, w1)
  w1 <- w1 / sum(w0, w1)
  plot(rbind(t(X0*w0), t(X1*w1)), asp = 1, pch = c("a", "b", "c", "t"), xlab = "x", ylab = "y")
  polygon(t(X0*w0), lty = 3)
})


# now with random data and more
dat <- simulate_data(N_covar = 2)
X0 <- dat$X0
X1 <- dat$X1

plot(rbind(t(X0), t(X1)), asp = 1, pch = c(rep("d", 50), "t"), xlab = "x", ylab = "y")

replicate(8, {
  w0 <- rexp(50)
  w1 <- rexp(1)
  plot(rbind(t(X0*w0), t(X1*w1)), asp = 1, pch = c(rep("d", 50), "t"), xlab = "x", ylab = "y")
})
