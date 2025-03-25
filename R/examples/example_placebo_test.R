set.seed(45)

# simulate data with an effect of 0.8 SD
dat <- simulate_data_synth(treatment_effect = .8)

# fit a model
fit <- pensynth(dat$X1, dat$X0, lambda = 1e-5)

# Perform placebo test
test <- placebo_test(fit, dat$Y1, dat$Y0)
plot(test)
abline(h = .8, lty = 2)
legend("bottomright", lty = 2, legend = "true effect")

# compute a pseudo p-value based on ATE in
# the post-intervention time period
ref_dist <- stats::ecdf(test$ATE0)
1 - ref_dist(test$ATE1)
