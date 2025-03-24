#' Simulate synthetic control data
#'
#' This function simulates a basic form of synthetic control
#' data, mainly for testing purposes.
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
    ar1_outcome = 0.8
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


#' Generate data from normal distribution with autocorrelation parameter
#'
#' @param n number of observations.
#' @param mean marginal mean
#' @param sd marginal standard deviation
#' @param phi autocorrelation parameter (-1 < phi < 1)
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
  if (abs(phi) >= 1)
    cli::cli_abort("AR parameter ({phi}) should be between -1 and 1.")
  if (phi == 0) return(rnorm(n, mean, sd))

  # create vector, pick first value, and fill remaining ones
  x <- numeric(n)
  x[1] <- rnorm(1, mean, sd)
  for (i in 2:n)
    x[i] <- rnorm(1, mean + phi*(x[i-1]-mean), sqrt(sd^2 - sd^2 * phi^2))
  return(x)
}

#' Generate data from independent multivariate normal
#' distribution with autocorrelation parameter
#'
#' @param n number of observations.
#' @param mean marginal mean
#' @param sd marginal standard deviation
#' @param phi autocorrelation parameter (-1 < phi < 1)
#'
#' @importFrom stats rnorm
#'
#' @return matrix of numeric values
#'
#' @details
#' Note that, unlike [stats::rnorm()], this function is
#' not vectorized over mean, sd, or phi.
#'
#' @keywords internal
mvrarnorm <- function(n, p, mean = 0, sd = 1, phi = 0) {
  # argument checks
  stopifnot(abs(phi) < 1)
  if (phi == 0) return(matrix(rnorm(n*p, mean, sd), n))

  # create vector, pick first value, and fill remaining ones
  x <- matrix(0.0, nrow = n, ncol = p)
  x[1,] <- rnorm(p, mean, sd)
  for (i in 2:n)
    x[i,] <- rnorm(p, mean + phi*(x[i-1,]-mean), sqrt(sd^2 - sd^2 * phi^2))
  return(x)
}

#' @export
simulate_data_factor <- function(
    N_donor = 50,
    N_nonzero = 4,
    N_covar = 5,
    N_pre = 12,
    N_post = 6,
    N_factors = 3,
    treatment_effect = 1,
    sd_factors = 1,
    ar1_factors = 0.8,
    sd_loadings = 1,
    sd_errors = 1,
    covar_means = FALSE
  ) {
  N_timepoints <- N_pre + N_post

  N_units <- N_donor + 1

  w  <- runif(N_donor)
  if (N_nonzero < N_donor) w[(N_nonzero + 1):N_donor] <- 0
  w  <- w / sum(w)

  # the factors, vary over time but not over units
  A <- mvrarnorm(N_timepoints, N_factors, sd = sd_factors, phi = ar1_factors)

  # the loadings for the outcome timeseries, vary over units but not over time
  L <- matrix(rnorm(N_donor*N_factors, sd = sd_loadings), N_factors)
  L <- cbind(L %*% w, L)

  # the residuals for outcome timeseries, vary over time and units
  E <- matrix(rnorm(N_timepoints*N_units, sd = sd_errors), N_timepoints)

  ZY <- A%*%L + E

  ZY[(N_pre + 1):N_timepoints, 1] <- ZY[(N_pre + 1):N_timepoints, 1] + treatment_effect

  # Now the covars
  for (i in 1:N_covar) {
    LCi <- matrix(rnorm(N_donor*N_factors, sd = sd_loadings), N_factors)
    LCi <- cbind(LCi %*% w, LCi)
    ECi <- matrix(rnorm(N_pre*N_units, sd = sd_errors), N_pre)
    Ci <- (A[1:N_pre,]%*%LCi + ECi)
    if (covar_means) {
      X <- if (i == 1) colMeans(Ci) else rbind(X, colMeans(Ci))
    } else {
      X <- if (i == 1) Ci else rbind(X, Ci)
    }
  }

  # Return matrices
  list(
    w = w,
    X0 = X[,-1], X1 = X[,1,drop=FALSE],
    Z0 = ZY[1:N_pre, -1], Z1 = ZY[1:N_pre,1, drop=FALSE],
    Y0 = ZY[(N_pre + 1):N_timepoints, -1], Y1 = ZY[(N_pre + 1):N_timepoints,1, drop = FALSE]
  )
}
