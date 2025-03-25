set.seed(45)
dat <- simulate_data_synth()
res <- with(dat, cv_pensynth(X1, X0, Z1, Z0))
plot(res)
