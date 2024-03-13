#' Penalized synthetic control estimator
#'
#' For a given set of variable weights (v) this function estimates
#' the unit weights for a synthetic control with penalization
#' according to Abadie & L'Hour (2021). This function deals with only a
#' single treated unit.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param v `N_covars vector` of variable weights
#' @param lambda `numeric` penalization parameter
#' @param opt_pars `clarabel` settings using [clarabel::clarabel_control()]
#' @param standardize `boolean` whether to standardize the input matrices (default TRUE)
#'
#' @details This routine uses the same notation of the original [Synth::synth()] implementation
#' but uses a different, faster quadratic program solver (namely, [osqp::osqp()]). Additionally, it
#' implements the penalization procedure described in Abadie & L'Hour (2021), such that the loss
#' function is as in equation 5 of that paper (but for a single treated unit).
#'
#' Variable weights are not optimized by this function, meaning they need to be pre-specified.
#' This is by design.
#'
#' The original synthetic control method can be recovered by setting lambda = 0. For determining
#' lambda based on data, see [cv_pensynth()].
#'
#' @references Abadie, A., & L’Hour, J. (2021).
#' A penalized synthetic control estimator for disaggregated data.
#' _Journal of the American Statistical Association, 116_(536), 1817-1834.
#'
#' @return A list with two values: `w`, the estimated weights; and
#' `solution`, the result of the optimization.
#'
#' @importFrom utils capture.output
#'
#' @examples
#' # generate some data
#' X0 <- matrix(
#'   c(1, 1.3,
#'     0.5, 1.8,
#'     1.1, 2.4,
#'     1.8, 1.8,
#'     1.3, 1.8), 2)
#' X1 <- matrix(c(0.8, 1.65), 2)
#' v <- rep(1, 2)
#'
#' # run classic synthetic control (no penalization)
#' res <- pensynth(X1, X0, v)
#' plot(t(X0))
#' points(t(X1), pch = 2)
#' points(t(X0%*%res$w), pch = 3)
#'
#' # run synthetic control with penalty
#' res <- pensynth(X1, X0, v, lambda = 0.5)
#' points(t(X0 %*% res$w), pch = 4)
#'
#' @seealso [cv_pensynth()] [Synth::synth()]
#'
#' @export
pensynth <- function(X1, X0, v, lambda = 0, opt_pars = clarabel::clarabel_control(), standardize = TRUE) {
  if (standardize) {
    st <- standardize_X(X1, X0)
    X0 <- st$X0
    X1 <- st$X1
  }
  N_donors <- ncol(X0)
  X0v <- X0*sqrt(v)
  X1v <- X1*sqrt(v)

  # components for quadratic program
  # see https://github.com/jeremylhour/pensynth/blob/master/functions/wsoll1.R
  X0VX0 <- crossprod(X0v)
  X1VX0 <- crossprod(X1v, X0v)
  Delta <- apply(X0v - c(X1v), 2, crossprod)

  # Constraint matrices
  Amat <- rbind(
    rep(1, N_donors), # Sum to 1 constraint
    -diag(N_donors) # Individ. weights gte 0 constraint
  )
  B <- c(
    1, # Sum to 1 constraint
    rep(0, N_donors) # Individ. weights gte 0 constraint
  )

# Run the quadratic program solver
  result <- clarabel::clarabel(
    P = X0VX0,
    q = -X1VX0 + lambda*Delta,
    A = Amat,
    b = B,
    cones = list(
      z = 1L, # There are 1 equalities
      l = N_donors # There are N_donors * 2 inequalities
    ),
    control = opt_pars
  )

  # clarabel only returns a numeric status code, so we'll add a
  # human-readable status column here (plus a description)
  result$status_description <- clarabel::solver_status_descriptions()[result$status][[1]]
  result$status <- names(clarabel::solver_status_descriptions()[result$status])


  return(list(w = result$x, solution = result))
}
