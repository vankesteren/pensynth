# Example of Bayesian bootstrap
set.seed(45)
dat <- simulate_data(N_donor = 50)
fit <- pensynth(dat$X1, dat$X0, lambda = .44)
bb  <- bootstrap(fit, n_iter = 500, X1 = dat$X1, X0 = dat$X0)

# Compute uncertainty interval
bb_pred <- predict(bb, rbind(dat$Z0, dat$Y0))
est     <- predict(fit, rbind(dat$Z0, dat$Y0))
upr95   <- apply(bb_pred, 1, quantile, 0.975)
lwr95   <- apply(bb_pred, 1, quantile, 0.025)

# create plot
plot(c(dat$Z1, dat$Y1), type = "l", ylim = range(c(dat$Z1, dat$Y1, upr95, lwr95)))
lines(est, col = "#343488")
polygon(c(1:18, 18:1), c(lwr95, rev(upr95)), col = "#34348855", border = NA)
abline(v = nrow(dat$Z1) + .5, lty = 2)


# plot posterior CI of ATE
ATE <- apply(predict(bb, dat$Y0), 2, function(Yhat) mean(dat$Y1 - Yhat))
plot(density(ATE), main = "ATE pseudo-posterior")
polygon(density(ATE), col = "#34348855", border = NA)
abline(v = mean(dat$Y1 - predict(fit, dat$Y0)), lty = 2)



# Verifying Bayesian bootstrap
library(pbapply)
cl <- parallel::makeCluster(10)
parallel::clusterEvalQ(cl, devtools::load_all())

ESIZE <- 0
NBOOT <- 50
parallel::clusterExport(cl, c("ESIZE", "NBOOT"))


# coverage test
res <- matrix(nrow = 100, ncol = 2)

sim_func <- function(i) {
  # gen data
  dli <- simulate_data(treatment_effect = ESIZE, N_donor = 50)
  fit <- pensynth(dli$X1, dli$X0, lambda = 1e-5)

  # Donor permutation test
  pt <- placebo_test(fit, Y1 = dli$Y1, Y0 = dli$Y0)
  res_pt <- ecdf(pt$ATE0)(pt$ATE1) > 0.975 | ecdf(pt$ATE0)(pt$ATE1) < 0.025

  # Bayesian bootstrap for ATE
  res_bb <- bootstrap(fit, n_iter = NBOOT, X1 = dli$X1, X0 = dli$X0)
  prd_bb <- predict(res_bb, dli$Y0)
  eff_bb <- sapply(1:NBOOT, \(i) dli$Y1 - prd_bb[,i])
  dst_bb <- ecdf(colMeans(eff_bb))
  lwr_bb <- quantile(dst_bb, probs = 0.025)
  upr_bb <- quantile(dst_bb, probs = 0.975)
  res_bb <- sign(lwr_bb * upr_bb) == 1
  return(c(res_pt, res_bb))
}

pbsapply(1:10, sim_func)




var_tot <- var(c(dli$Z1 - predict(fit, dli$Z0))) + apply(eff_bb, 1, var)

lwr <- mean(dli$Y1 - predict(fit, dli$Y0) - 1.96 * sqrt(var_tot))
upr <- mean(dli$Y1 - predict(fit, dli$Y0) + 1.96 * sqrt(var_tot))

eff_bb


# some SCPI comparison
compute_bayesian_bootstrap <- function(dat) {
  res_bb <- sapply(1:NBOOT, \(i) pensynth(dat$X1, dat$X0, lambda = 1e-5, boot_weight = TRUE)$w)
  eff_bb <- sapply(1:NBOOT, \(i) dat$Y1 - dat$Y0 %*% res_bb[,i])
  dst_bb <- ecdf(colMeans(eff_bb))
  lwr_bb <- quantile(dst_bb, probs = 0.025)
  upr_bb <- quantile(dst_bb, probs = 0.975)
  cvg_bb <- lwr_bb < ESIZE && upr_bb > ESIZE
  len_bb <- upr_bb - lwr_bb
  return(c(cvg_bb, len_bb))
}

scpi_format <- function(dat) {
  scpi_specs <- list(
    J = ncol(dat$X0),
    K = rep(0, nrow(dat$X0)),
    KM = 0,
    M = nrow(dat$X0),
    I = 1,
    cointegrated.data = FALSE,
    period.pre = 1:nrow(dat$Z0),
    period.post = (1:nrow(dat$Y0)) + nrow(dat$Z0),
    T0.features = rep(1, nrow(dat$X0)),
    T1.outcome = nrow(dat$Y0),
    outcome.var = "outcome",
    features = paste0("X", 1:nrow(dat$X0)),
    constant = FALSE,
    out.in.features = FALSE,
    treated.units = "treated unit",
    donors.units = paste("donor", 1:ncol(dat$X0)),
    effect = "unit-time",
    sparse.matrices = FALSE,
    units.est = "treated unit"
  )
  attr(scpi_specs$K, "names") <- attr(scpi_specs$T0.features, "names") <- scpi_specs$features
  # do dimnames
  colnames(dat$X0) <- colnames(dat$Y0) <- colnames(dat$Z0) <-
    paste(scpi_specs$treated.units, scpi_specs$donors.units, sep = ".")
  rownames(dat$X1) <- rownames(dat$X0) <- paste(scpi_specs$treated.units, scpi_specs$features, 1, sep = ".")
  rownames(dat$Y0) <- rownames(dat$Y1) <- paste(scpi_specs$treated.units, scpi_specs$period.post, sep = ".")
  rownames(dat$Z0) <- rownames(dat$Z1) <- paste(scpi_specs$treated.units, scpi_specs$period.pre, sep = ".")
  structure(list(
    A = dat$X1,
    B = dat$X0,
    C = NULL,
    P = dat$Y0,
    Y.pre = dat$Z1,
    Y.post = dat$Y1,
    Y.donors = dat$Z0,
    specs = scpi_specs
  ), class = "scdata")
}

scpi_dat <- scpi_format(dat)
res <- scpi(scpi_dat, w.constr = list(name = "simplex", Q = 1), e.order = 0, u.order = 0, rho = .44)
res$inference.results$CI.all.qreg
scplot(res)$plot_out_qr

# scpi plot but using bb
resid_pre  <- apply(predict(bb, dat$Z0), 2, \(x) x - dat$Z1)
# https://stats.stackexchange.com/a/501275/116878
fuzzed_pred <- predict(bb, dat$Y0) + sample(resid_pre, size = length(bb$boot_fit)*nrow(dat$Y0), replace = TRUE)

apply(fuzzed_pred, 1, quantile, prob = 0.025)
apply(fuzzed_pred, 1, quantile, prob = 0.975)
pre_pred <- predict(bb, dat$Z0)
data.frame(
  estim = predict(fit, rbind(dat$Z0, dat$Y0)),
  lower = apply(rbind(pre_pred, fuzzed_pred), 1, quantile, probs = 0.025),
  upper = apply(rbind(pre_pred, fuzzed_pred), 1, quantile, probs = 0.975)
) |> mutate(width = upper - lower) |>
  mutate(time = 1:length(c(dat$Z1, dat$Y1)), pre = time <= length(dat$Z1)) |>
  ggplot(aes(x = time, y = estim, ymin = lower, ymax = upper)) +
  geom_vline(xintercept = nrow(dat$Z0) + 0.5, linetype = 2) +
  geom_line(linetype = 2) +
  geom_line(aes(y = rowMeans(predict(bb, rbind(dat$Z0, dat$Y0))))) +
  geom_line(aes(colour = pre), linetype = 2) +
  geom_errorbar(aes(colour = pre), width = 0.2) +
  geom_point(aes(colour = pre), size = 2.5) +
  geom_line(aes(y = c(dat$Z1, dat$Y1))) +
  geom_point(aes(y = c(dat$Z1, dat$Y1)), pch = 21, size = 2.5, stroke = .8) +
  scale_colour_manual(values = c("TRUE" = "seagreen", "FALSE" = "blue"), guide = "none") +
  theme_linedraw() +
  labs(
    x = "timepoint",
    y = "outcome",
    title = "Bayesian bootstrap pseudo-prediction intervals",
    subtitle = paste(length(bb$boot_fit), "bootstrap replications")
  )


# test boot coverage vs scpi coverage
# seems much better????
sim_func <- function(effect = 0) {
  # simulate some data
  dat <- simulate_data(treatment_effect = effect)

  # scpi coverage
  dat_sc <- scpi_format(dat)

  scpi_res <- suppressWarnings(
    scpi(
      data = dat_sc,
      w.constr = list(name = "simplex", Q = 1),
      e.order = 0,
      u.order = 0,
      rho = 1e-5,
      verbose = FALSE
    )
  )
  scpi_eff <- colMeans(apply(scpi_res$inference.results$CI.all.qreg, 2, \(Yhat) dat$Y1 - Yhat))[1:2]

  # compute boot solution
  fit <- pensynth(dat$X1, dat$X0, lambda = 1e-5)
  boo <- bootstrap(fit, n_iter = 200L, X1 = dat$X1, X0 = dat$X0)

  # compute residuals
  pre <- predict(boo, dat$Z0)
  rsd <- apply(pre, 2, function(Zhat) dat$Z1 - Zhat)

  # create predictions w/ boot residuals
  prd <- predict(boo, dat$Y0) + sample(rsd, nrow(dat$Y0) * length(boo), replace = TRUE)
  eff <- apply(prd, 2, function(Yhat) dat$Y1 - Yhat)
  ate <- colMeans(eff)
  lwr <- quantile(ate, 0.025)
  upr <- quantile(ate, 0.975)

  c("scpi_lo" = scpi_eff[2], "scpi_hi" = scpi_eff[1], "boot_lo" = lwr, "boot_hi" = upr)
}

clus <- makeCluster(12)
clusterEvalQ(clus, devtools::load_all())
clusterEvalQ(clus, library(scpi))
clusterExport(clus, c("sim_func", "scpi_format"))
res <- pbsapply(1:100, \(i) sim_func(1.56), cl = clus)

library(tidyverse)
t(res) |>
  as_tibble() |>
  pivot_longer(everything(), names_to = c("method", ".value"), names_pattern = "(.*)_([hilo]{2})") |>
  mutate(width = hi - lo, cvg = hi > 1.56 & lo < 1.56) |>
  summarize(across(everything(), mean), .by = method)

rowMeans(res)

library(pbapply)
res <- pbsapply(1:1000, function(i) {
  dati <- simulate_data()
  fiti <- pensynth(dati$X1, dati$X0)
  booi <- bootstrap(fiti, n_iter = 200L, X1 = dati$X1, X0 = dati$X0)
  prdi <- predict(booi, dati$Y0)
  mean(prdi)
}, cl = clus)
hist(res)


boot_intervals <- function(boo, Z1, Z0, Y1, Y0, level = 0.95) {
  cutoffs <- c(.5 - level / 2, .5 + level / 2)

  # First, predict pre-intervention outcome
  pred_Z1 <- predict(boo, Z0)
  smp_int_Z1 <- t(apply(pred_Z1, 1, quantile, probs = cutoffs))


  # compute bootstrapped prediction residuals
  residuals <- c(apply(pred_Z1, 2, function(Z1hat) Z0 - Z1hat))

  # prediction intervals are sampling uncertainty + prediction uncertainty
  fuzz_Z1 <- pred_Z1 + sample(residuals, length(boo) * nrow(Z0), replace = TRUE)
  prd_int_Z1  <- t(apply(fuzz_Z1, 1, quantile, probs = cutoffs))

  pred_Y1 <- predict(boo, Y0)
  smp_int_Y1 <- t(apply(pred_Y1, 1, quantile, probs = cutoffs))
  fuzz_Y1 <- pred_Y1 + sample(residuals, length(boo) * nrow(Y0), replace = TRUE)
  prd_int_Y1 <- t(apply(fuzz_Y1, 1, quantile, probs = cutoffs))

  # compute effect intervals as well
  pred_Eff <- apply(pred_Y1, 2, function(Y1hat) Y1 - Y1hat)
  smp_int_Eff <- t(apply(pred_Eff, 1, quantile, probs = cutoffs))
  fuzz_Eff <- apply(fuzz_Y1, 2, function(Y1hat) Y1 - Y1hat)
  prd_int_Eff <- t(apply(fuzz_Eff, 1, quantile, probs = cutoffs))

  # add colnames
  colnames(smp_int_Z1)  <- c("lower", "upper")
  colnames(prd_int_Z1)  <- c("lower", "upper")
  colnames(smp_int_Y1)  <- c("lower", "upper")
  colnames(prd_int_Y1)  <- c("lower", "upper")
  colnames(smp_int_Eff) <- c("lower", "upper")
  colnames(prd_int_Eff) <- c("lower", "upper")

  return(list(
    samp_interval_pre    = smp_int_Z1,
    pred_interval_pre    = prd_int_Z1,
    samp_interval_post   = smp_int_Y1,
    pred_interval_post   = prd_int_Y1,
    samp_interval_effect = smp_int_Eff,
    pred_interval_effect = prd_int_Eff
  ))
}

dat <- simulate_data(treatment_effect = 3)
fit <- cv_pensynth(dat$X1, dat$X0, dat$Z1[1:6*2-1,], dat$Z0[1:6*2-1,])
boo <- bootstrap(fit, n_iter = 200L, dat$X1, dat$X0, dat$Z1[1:6*2-1,], dat$Z0[1:6*2-1,])
int_pre  <- boot_intervals(boo, dat$Z1, dat$Z0, dat$Y1, dat$Y0, level = 0.95)
int_post <- boot_intervals(boo, dat$Z1[1:6*2,], dat$Z0[1:6*2,], dat$Y1, dat$Y0, level = 0.95)

library(tidyverse)
as_tibble(rbind(int_pre$samp_interval_pre, int_post$pred_interval_post)) |>
  mutate(
    est = c(predict(fit, dat$Z0), predict(fit, dat$Y0)),
    time = 1:n(),
    post = time > nrow(dat$Z0)
  ) |>
  ggplot(aes(x = time, y = est, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2, fill = "seagreen") +
  geom_line(linetype = 2) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(y = c(dat$Z1, dat$Y1))) +
  geom_point(aes(y = c(dat$Z1, dat$Y1)), size = 2) +
  geom_vline(xintercept = nrow(dat$Z0) + 0.5, linetype = 2) +
  theme_linedraw() +
  #scale_alpha_manual(values = c(0.1, 0.5), guide = "none") +
  labs(
    x = "Timepoint",
    y = "Outcome",
    title = "Prediction"
  )
