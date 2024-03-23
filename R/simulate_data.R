#' Simulate synthetic control data
#'
#' This function simulates a basic form of synthetic control
#' data, mainly for testing purposes. This
#'
#' @param N_donor number of donors
#' @param N_covar number of covariates
#' @param N_pre number of pre-intervention timepoints
#' @param N_post number of post-intervention timepoints
#' @param N_nonzero number of true nonzero weights
#' @param treatment_effect the size of the true treatment effect
#' @param sd_resid_X1 the residual standard deviation of X1
#' @param sd_resid_Z1 the residual standard deviation of Z1
#' @param sd_resid_Y1 the residual standard deviation of Y1
#'
#' @returns A list with the following elements
#' - w the true unit weights
#' - X0 the donor unit covariates
#' - X1 the treated unit covariates
#' - Z0 the donor unit pre-intervention outcomes
#' - Z1 the treated unit pre-intervention outcomes
#' - Y0 the donor unit post-intervention outcomes
#' - Y1 the treated unit post-intervention outcomes
#'
#' @details
#' Note that treatment effect can be a single number, but
#' it may also be a vector of length N_post, indicating
#' the effect size at each post-intervention measurement
#' occasion.
#'
#' @export
simulate_data <- function(
    N_donor = 50,
    N_covar = 5,
    N_pre = 12,
    N_post = 6,
    N_nonzero = 4,
    treatment_effect = 1,
    sd_resid_X1 = 0.1,
    sd_resid_Z1 = 0.1,
    sd_resid_Y1 = 0.1
  ) {
  w  <- runif(N_donor)
  if (N_nonzero < N_donor) w[(N_nonzero + 1):N_donor] <- 0
  w  <- w / sum(w)
  X0 <- matrix(rnorm(N_covar * N_donor), N_covar)
  X1 <- X0 %*% w + rnorm(N_covar, sd = sd_resid_X1)
  Z0 <- matrix(rnorm(N_pre * N_donor), N_pre)
  Z1 <- Z0 %*% w + rnorm(N_pre, sd = sd_resid_Z1)
  Y0 <- matrix(rnorm(N_post * N_donor), N_post)
  Y1 <- Y0 %*% w + treatment_effect + rnorm(N_post, sd = sd_resid_Y1)
  list(w = w, X0 = X0, X1 = X1, Z0 = Z0, Z1 = Z1, Y0 = Y0, Y1 = Y1)
}
