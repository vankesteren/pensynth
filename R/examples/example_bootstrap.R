# Verifying Bayesian bootstrap
library(pbapply)
devtools::load_all()
cl <- parallel::makeCluster(10)
parallel::clusterEvalQ(cl, devtools::load_all())

ESIZE <- 0.5
NBOOT <- 2000

dat <- simulate_data(treatment_effect = ESIZE, N_donor = 50)
parallel::clusterExport(cl, "dat")
res_bb <- pbsapply(1:NBOOT, \(i) {
  cv_pensynth(dat$X1, dat$X0, dat$Z1, dat$Z0, boot_weight = TRUE)$w_opt
}, cl = cl)

Y1_bb <- pbsapply(1:NBOOT, \(i) dat$Y0 %*% res_bb[,i])
Z1_bb <- pbsapply(1:NBOOT, \(i) dat$Z0 %*% res_bb[,i])

par(mfrow = c(2, 1))
# prediction plot
plot(NA, ylim = c(-3, 3), xlim = c(1, 18), main = "Prediction plot", xlab = "outcome")
for (i in 1:NBOOT) lines(c(Z1_bb[,i], Y1_bb[,i]), col = "#6666AA11")
lines(c(rowMeans(Z1_bb), rowMeans(Y1_bb)), lty = 2, lwd = 2, col = "#000088")
lines(c(dat$Z1, dat$Y1), lwd = 2)
abline(v = 12.5, lty = 3)

# effect plot
eff_bb <- sapply(1:NBOOT, \(i) dat$Y1 - dat$Y0 %*% res_bb[,i])
hist(colMeans(eff_bb), breaks = 70, main = "ATE pseudo-posterior", freq = FALSE)
abline(v = mean(eff_bb), lty = 3, lwd = 2)
abline(v = ESIZE, lwd = 2)
par(mfrow = c(1, 1))

# bootstrap p-value
ref_dist <- ecdf(colMeans(eff_bb))
if (ref_dist(0) > 0.5) 1 - ref_dist(0) else ref_dist(0)


# coverage test
parallel::clusterExport(cl, c("ESIZE", "NBOOT"))

res <- pbsapply(1:1000, function(yoyoyo) {
  # gen data
  dat <- simulate_data(treatment_effect = ESIZE, N_donor = 20)

  # Bayesian bootstrap
  res_bb <- sapply(1:NBOOT, \(i) pensynth(dat$X1, dat$X0, lambda = 1e-5, boot_weight = TRUE)$w)
  eff_bb <- sapply(1:NBOOT, \(i) dat$Y1 - dat$Y0 %*% res_bb[,i])
  dst_bb <- ecdf(colMeans(eff_bb))
  lwr_bb <- quantile(dst_bb, probs = 0.025)
  upr_bb <- quantile(dst_bb, probs = 0.975)
  cvg_bb <- lwr_bb < ESIZE && upr_bb > ESIZE
  len_bb <- upr_bb - lwr_bb


  # default bootstrap
  res_bs <- sapply(1:NBOOT, \(i) pensynth(dat$X1, dat$X0[,sample(20, replace = TRUE)], lambda = 1e-5, boot_weight = TRUE)$w)
  eff_bs <- sapply(1:NBOOT, \(i) dat$Y1 - dat$Y0 %*% res_bs[,i])
  dst_bs <- ecdf(colMeans(eff_bs))
  lwr_bs <- quantile(dst_bs, probs = 0.025)
  upr_bs <- quantile(dst_bs, probs = 0.975)
  cvg_bs <- lwr_bs < ESIZE && upr_bs > ESIZE
  len_bs <- upr_bs - lwr_bs

  return(c(cvg_bb, len_bb, cvg_bs, len_bs))
}, cl = cl)

rowMeans(res)

parallel::stopCluster(cl)


