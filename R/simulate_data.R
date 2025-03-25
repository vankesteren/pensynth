#' Simulate data according to synthetic control model
#'
#' This function simulates a basic form of synthetic control
#' data, mainly for testing purposes.
#'
#' @param N_donor number of donors
#' @param N_treated number of treated units
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
#' - W the true unit weights
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
#' the effect size at each post-intervention measurement.
#' occasion. It may also be a matrix of size N_post by
#' N_treated.
#'
#'
#' @seealso [pensynth()], [cv_pensynth()], [placebo_test()], [simulate_data_factor()]
#'
#' @export
simulate_data <- function(N_donor = 50,
                          N_treated = 1,
                          N_covar = 5,
                          N_pre = 12,
                          N_post = 6,
                          N_nonzero = 4,
                          treatment_effect = 1,
                          sd_resid_X = 0.1,
                          sd_resid_ZY = 0.1,
                          ar1_outcome = 0.8) {
  # Total units
  N_units <- N_donor + N_treated

  # Weights
  W <- matrix(runif(N_donor * N_treated), N_donor)
  if (N_nonzero < N_donor)
    for (n in 1:N_treated) {
      id_zero <- sample(N_donor, N_donor - N_nonzero, replace = FALSE)
      W[id_zero, n] <- 0
    } else {
      cli::cli_alert_warning("N_nonzero ({N_nonzero}) larger than N_donor ({N_donor})")
    }
  W <- apply(W, 2, function(w) w / sum(w))

  # Covariates
  X0 <- matrix(rnorm(N_covar * N_donor), N_covar)
  X1 <- X0 %*% W + matrix(rnorm(N_covar * N_treated, sd = sd_resid_X), N_covar)

  # Outcome
  N_tot <- N_pre + N_post

  # Pre- and post outcome for donor units
  ZY0 <- mvrarnorm(N_tot, N_donor, phi = ar1_outcome)
  Z0  <- ZY0[1:N_pre, ]
  Y0  <- ZY0[(N_pre + 1):N_tot, ]

  # Pre- and post outcome for treated units
  RZY <- matrix(rnorm(N_tot * N_treated, sd = sd_resid_ZY), N_tot)
  RZ  <- RZY[1:N_pre, ]
  RY  <- RZY[(N_pre + 1):N_tot, ]

  Z1 <- Z0 %*% W + RZ
  Y1 <- Y0 %*% W + treatment_effect + RY

  # Return matrices
  list(
    W = W,
    X0 = X0,
    X1 = X1,
    Z0 = Z0,
    Z1 = Z1,
    Y0 = Y0,
    Y1 = Y1
  )
}

#' Simulate data according to factor model
#'
#' This function simulates data according to a latent factor model:
#' 1. Simulate time-varying latent factors, which are the same for all units
#' 2. Simulate time-invariant factor loadings, separately for each donor unit
#' 3. Create sparse unit weights for each treated unit
#' 3. Compute the loadings for the treated units as donor-unit loadings * weights
#' 4. Simulate observed outcome time-series as factors * loadings + error
#' 5. Do the same for each covariate, holding loadings equal. Average across pre-intervention timepoints.
#'
#' @param N_donor number of donors
#' @param N_treated number of treated units
#' @param N_nonzero number of true nonzero weights
#' @param N_covar number of covariates
#' @param N_pre number of pre-intervention timepoints
#' @param N_post number of post-intervention timepoints
#' @param N_factors number of latent factors to simulate
#' @param treatment_effect the size of the true treatment effect
#' @param sd_factors the standard deviation of the (unit-invariant, time-varying) factors
#' @param ar1_factors autoregressive effect of the factors
#' @param sd_loadings the standard deviation of the (time-invariant) factor loadings
#' @param sd_errors the standard deviation of the independent errors
#' @param covar_means whether to average the covariates across the pre-intervention times (experimental)
#'
#' @returns A list with the following elements
#' - W the true unit weights
#' - X0 the donor unit covariates
#' - X1 the treated unit covariates
#' - Z0 the donor unit pre-intervention outcomes
#' - Z1 the treated unit pre-intervention outcomes
#' - Y0 the donor unit post-intervention outcomes
#' - Y1 the treated unit post-intervention outcomes
#'
#' @importFrom stats rnorm runif
#'
#' @example R/examples/example_simulate_factor.R
#'
#' @details
#' Note that treatment effect can be a single number, but
#' it may also be a vector of length N_post, indicating
#' the effect size at each post-intervention measurement.
#' occasion. It may also be a matrix of size N_post by
#' N_treated.
#'
#' Standard values of sd_factors, sd_loadings, and sd_errors
#' have been chosen such that the observed variables have
#' expected variance of 1.
#'
#'
#' @seealso [pensynth()], [cv_pensynth()], [placebo_test()], [simulate_data()]
#'
#' @export
simulate_data_factor <- function(N_donor = 50,
                                 N_treated = 1,
                                 N_nonzero = 4,
                                 N_covar = 5,
                                 N_pre = 12,
                                 N_post = 6,
                                 N_factors = 3,
                                 treatment_effect = 1,
                                 sd_factors = sqrt(2) / 2,
                                 ar1_factors = 0.8,
                                 sd_loadings = sqrt(2) / 2,
                                 sd_errors = 0.5,
                                 covar_means = TRUE) {
  N_timepoints <- N_pre + N_post

  N_units <- N_donor + N_treated

  W <- matrix(runif(N_donor * N_treated), N_donor)
  if (N_nonzero < N_donor)
    for (n in 1:N_treated) {
      id_zero <- sample(N_donor, N_donor - N_nonzero, replace = FALSE)
      W[id_zero, n] <- 0
    } else {
      cli::cli_alert_warning("N_nonzero ({N_nonzero}) larger than N_donor ({N_donor})")
    }
  W <- apply(W, 2, function(w)
    w / sum(w))

  # the factors, vary over time but not over units
  A <- mvrarnorm(N_timepoints, N_factors, sd = sd_factors, phi = ar1_factors)

  # the loadings for the outcome timeseries, vary over units but not over time
  L <- matrix(rnorm(N_donor * N_factors, sd = sd_loadings), N_factors)

  # add treated unit loadings
  L <- cbind(L %*% W, L)

  # the residuals for outcome timeseries, vary over time and units
  E <- matrix(rnorm(N_timepoints * N_units, sd = sd_errors), N_timepoints)

  # compute Z (pre-intervention timeseries) and Y (post-intervention timeseries)
  # in absence of treatment effect (we will add this at the end)
  ZY <- A %*% L + E

  # Now the covars
  X <- vector("list", N_covar)
  for (i in 1:N_covar) {
    # same exact process as outcome timeseries
    LCi <- matrix(rnorm(N_donor * N_factors, sd = sd_loadings), N_factors)
    LCi <- cbind(LCi %*% W, LCi)
    ECi <- matrix(rnorm(N_pre * N_units, sd = sd_errors), N_pre)
    Ci <- (A[1:N_pre, ] %*% LCi + ECi)
    if (covar_means) {
      X <- if (i == 1)
        colMeans(Ci)
      else
        rbind(X, colMeans(Ci))
    } else {
      X <- if (i == 1)
        Ci
      else
        rbind(X, Ci)
    }
  }

  # Return matrices
  list(
    W  = W,
    X0 = unname(X[, (N_treated + 1):N_units, drop = FALSE]),
    X1 = unname(X[, 1:N_treated, drop = FALSE]),
    Z0 = ZY[1:N_pre, (N_treated + 1):N_units, drop = FALSE],
    Z1 = ZY[1:N_pre, 1:N_treated, drop = FALSE],
    Y0 = ZY[(N_pre + 1):N_timepoints, (N_treated + 1):N_units],
    Y1 = ZY[(N_pre + 1):N_timepoints, 1:N_treated, drop = FALSE] + treatment_effect
  )
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
rarnorm <- function(n,
                    mean = 0,
                    sd = 1,
                    phi = 0) {
  # argument checks
  if (abs(phi) >= 1)
    cli::cli_abort("AR parameter ({phi}) should be between -1 and 1.")
  if (phi == 0)
    return(rnorm(n, mean, sd))

  # create vector, pick first value, and fill remaining ones
  x <- numeric(n)
  x[1] <- rnorm(1, mean, sd)
  for (i in 2:n)
    x[i] <- rnorm(1, mean + phi * (x[i - 1] - mean), sqrt(sd^2 - sd^2 * phi^2))
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
mvrarnorm <- function(n,
                      p,
                      mean = 0,
                      sd = 1,
                      phi = 0) {
  # argument checks
  stopifnot(abs(phi) < 1)
  if (phi == 0)
    return(matrix(rnorm(n * p, mean, sd), n))

  # create vector, pick first value, and fill remaining ones
  x <- matrix(0.0, nrow = n, ncol = p)
  x[1, ] <- rnorm(p, mean, sd)
  for (i in 2:n)
    x[i, ] <- rnorm(p, mean + phi * (x[i - 1, ] - mean), sqrt(sd^2 - sd^2 * phi^2))
  return(x)
}
