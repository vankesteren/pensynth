# Showcasing features of the pensynth package
set.seed(45)

# simulate data with an effect of 0.8 SD
dat <- simulate_data(treatment_effect = 0.8, N_donor = 150)
plot(c(dat$Z1, dat$Y1), type = "l", ylim = c(-2.5, 2.5),
     main = "Observed data")
abline(v = length(dat$Z1) + 0.5, lty = 3)

# is the treated unit in the convex hull of the donors?
in_convex_hull(dat$X1, dat$X0)

# fit a hold-out-validated penalized synthetic control
fit <- cv_pensynth(dat$X1, dat$X0, dat$Z1, dat$Z0)

# plot the solution path with hold-out MSE
plot(fit)

# compute and plot estimated synthetic control
Z1hat <- predict(fit, dat$Z0)
Y1hat <- predict(fit, dat$Y0)
plot(c(dat$Z1, dat$Y1), type = "l", ylim = c(-2.5, 2.5),
     main = "Synthetic control plot")
lines(c(Z1hat, Y1hat), lty = 2)
abline(v = length(Z1hat) + 0.5, lty = 3)
legend(
  "topleft",
  legend = c("Observed outcome", "Synthetic control", "Intervention"),
  lty = 1:3
)

# check if true effect of 0.8 is estimated well
mean(dat$Y1 - Y1hat) # 0.886, not bad!

# Perform placebo (permutation) test
test <- placebo_test(fit, dat$Y1, dat$Y0)
plot(test)
abline(h = .8, lty = 2)
legend("bottomright", lty = 2, legend = "true effect")

# compute a pseudo p-value based on ATE in
# the post-intervention time period
ref_dist <- stats::ecdf(test$ATE0)
1 - ref_dist(test$ATE1)
