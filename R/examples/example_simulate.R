# simulate data with an effect of 0.8 SD
dat <- simulate_data(treatment_effect = 0.8)

plot(
  NA,
  ylim = c(-3, 3),
  xlim = c(1, 18),
  main = "Simulated data",
  ylab = "Outcome value"
)
for (n in 1:ncol(dat$Z0))
  lines(1:18, c(dat$Z0[, n], dat$Y0[, n]), col = "grey")
lines(1:18, c(dat$Z1, dat$Y1))
lines(1:18, rbind(dat$Z0, dat$Y0) %*% dat$w, lty = 2)
abline(v = length(dat$Z1) + 0.5, lty = 3)
legend(
  x = "bottomleft",
  legend = c(
    "Donor units",
    "Treated unit",
    "True synth. control",
    "Intervention time"
  ),
  lty = c(1, 1, 2, 3),
  col = c("grey", "black", "black", "black")
)
