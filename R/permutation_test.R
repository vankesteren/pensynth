permutation_test <- function(object, Y1, Y0, plot = TRUE) {
  # original obs - pred
  est <- Y1 - predict(object, Y0)

  # create a list of weights
  N_donors <- length(object$w)
  ps_list <- lapply(
    1:N_donors,
    function(n) {
      new_X1 <- eval(getCall(object)$X0)[, n, drop = FALSE]
      new_X0 <- cbind(
        eval(getCall(object)$X1),
        eval(getCall(object)$X0)[, -n, drop = FALSE]
      )
      update(object = object, X1 = new_X1, X0 = new_X0)
    }
  )

  # create prediction for each iter
  null <- sapply(
    1:N_donors,
    function(n) {
      Y0[, n, drop = FALSE] - predict(ps_list[[n]], cbind(Y1, Y0[, -n, drop = FALSE]))
    }
  )

  # do plot if needed
  out <- list(E0 = null_effects, E1 = est_effect)
  if (plot) {
    plot(NA, ylim = range(c(null, est)), xlim = c(1, nrow(null)),
         ylab = "Treatment effect")
    for (n in 1:N_donors) lines(null[,n], col = "grey")
    lines(est)
    legend("topright", lty = c(1, 1), col = c("grey", "black"),
           legend = c("reference", "estimate"))
  }

  return(out)
}

# dat <- sim_dat(N_donor = 60, sdx = .5, sdz = .5, effect = 3/(1:6))
# res  <- pensynth(dat$X1, dat$X0, 1, 0.1)
# pres <- permutation_test(res, dat$Y1, dat$Y0)
# pval <- 1 - ecdf(apply(pres$E0, 2, crossprod))(crossprod(pres$E1))
# pval
