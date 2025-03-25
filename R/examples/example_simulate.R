# simulate data with an effect of 0.8 SD
dat <- simulate_data_synth(treatment_effect = 0.8)

plot(
  NA,
  ylim = c(-3, 3),
  xlim = c(1, 18),
  main = "Simulated data",
  ylab = "Outcome value",
  xlab = "Timepoint"
)
for (n in 1:ncol(dat$Z0))
  lines(1:18, c(dat$Z0[, n], dat$Y0[, n]), col = "grey")
lines(1:18, c(dat$Z1, dat$Y1), lwd = 2)
lines(1:18, rbind(dat$Z0, dat$Y0) %*% dat$W, lty = 2, lwd = 2)
abline(v = length(dat$Z1) + 0.5, lty = 3)
legend(
  x = "bottomleft",
  legend = c(
    "Donor units",
    "Treated unit",
    "Synth. control"
  ),
  lty = c(1, 1, 2),
  lwd = c(1, 2, 2),
  col = c("grey", "black", "black")
)
text(length(dat$Z1) + 0.5, -3, "Intervention\ntimepoint", pos = 4, font = 3)
