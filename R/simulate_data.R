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
#' @param sd_resid_X the residual standard deviation of X1
#' @param sd_resid_ZY the residual standard deviation of Z1 and Y1
#' @param ar1_outcome autoregressive effect of the outcome
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
#' @importFrom stats rnorm runif
#'
#' @example R/examples/example_simulate.R
#'
#' @details
#' Note that treatment effect can be a single number, but
#' it may also be a vector of length N_post, indicating
#' the effect size at each post-intervention measurement
#' occasion.
#'
#'
#'
#' @seealso [pensynth()], [cv_pensynth()], [placebo_test()]
#'
#' @export
simulate_data <- function(
    N_donor = 50,
    N_covar = 5,
    N_pre = 12,
    N_post = 6,
    N_nonzero = 4,
    treatment_effect = 1,
    sd_resid_X = 0.1,
    sd_resid_ZY = 0.1,
    ar1_outcome = 0
  ) {

  # Weights
  w  <- runif(N_donor)
  if (N_nonzero < N_donor) w[(N_nonzero + 1):N_donor] <- 0
  w  <- w / sum(w)

  # Covariates
  X0 <- matrix(rnorm(N_covar * N_donor), N_covar)
  X1 <- X0 %*% w + rnorm(N_covar, sd = sd_resid_X)

  # Outcome
  N_tot <- N_pre + N_post

  # Pre- and post outcome for donor units
  ZY0 <- matrix(rarnorm(N_tot * N_donor, phi = ar1_outcome), N_tot)
  Z0  <- ZY0[1:N_pre,]
  Y0  <- ZY0[(N_pre+1):N_tot,]

  # Pre- and post outcome for treated unit
  RZY <- matrix(rnorm(N_tot * 1, sd = sd_resid_ZY), N_tot)
  RZ  <- RZY[1:N_pre,]
  RY  <- RZY[(N_pre+1):N_tot,]

  Z1 <- Z0 %*% w + RZ
  Y1 <- Y0 %*% w + treatment_effect + RY

  # Return matrices
  list(w = w, X0 = X0, X1 = X1, Z0 = Z0, Z1 = Z1, Y0 = Y0, Y1 = Y1)
}


#' Generate data from normal distribution with AR1 parameter
#'
#' @param n number of observations.
#' @param mean marginal mean
#' @param sd marginal standard deviation
#' @param phi autoregressive parameter (-1 <= phi <= 1)
#'
#' @importFrom stats rnorm
#'
#' @return vector of numeric values
#'
#' @details
#' Note that, unlike [stats::rnorm()], this function is
#' not vectorized over mean, sd, or phi.
#'
#' @keywords internal
rarnorm <- function(n, mean = 0, sd = 1, phi = 0) {
  # argument checks
  stopifnot(phi >= -1 & phi <= 1)
  if (phi == 0) return(rnorm(n, mean, sd))

  # create vector, pick first value, and fill remaining ones
  x <- numeric(n)
  x[1] <- rnorm(1, mean, sd)
  for (i in 2:n)
    x[i] <- rnorm(1, mean + phi*(x[i-1]-mean), sqrt(sd^2 - sd^2 * phi^2))
  return(x)
}
