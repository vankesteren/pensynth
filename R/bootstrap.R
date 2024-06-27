#' Bayesian bootstrap inference
#'
#' Produce Bayesian bootstrap prediction intervals and ATE intervals
#' for pensynth objects. This is done through reweighing the columns
#' of the data (X0 and X1) using Dirichlet weights.
#'
#' @param object a fitted `pensynth` or `cvpensynth` object
#' @param n_iter `int` number of bootstrap iterations (default 1000)
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param Z1 `N_targets by 1 matrix` of treated unit hold-out outcome
#' @param Z0 `N_targets by N_donors matrix` of donor unit hold-out outcome
#'
#' @export
bootstrap <- function(object, n_iter = 200L, ...) {
  UseMethod("bootstrap")
}

#' @rdname bootstrap
#' @export
bootstrap.pensynth <- function(object, n_iter, X1, X0) {
  nweight <- ncol(X0) + ncol(X1)
  # create reweighed predictions
  bb_fit <- lapply(
    X = 1:n_iter,
    FUN = function(n) {
      # rerun original call but with boot weights
      cat("\rProgress: [", n, "/", n_iter, "]")
      update(object = object, X1 = X1, X0 = X0, iweights = rexp(nweight))
    }
  )
  cat("\n")
  return(structure(
    .Data = bb_fit,
    class = "pensynthboot"
  ))
}

#' @rdname bootstrap
#' @export
bootstrap.cvpensynth <- function(object, n_iter, X1, X0, Z1, Z0) {
  nweight <- ncol(X0) + ncol(X1)
  # create reweighed predictions
  bb_fit <- lapply(
    X = 1:n_iter,
    FUN = function(n) {
      # rerun original call but with boot weights
      cat("\rProgress: [", n, "/", n_iter, "]")
      update(object = object, X1 = X1, X0 = X0, Z1 = Z1, Z0 = Z0, iweights = rexp(nweight))
    }
  )
  cat("\n")
  return(structure(
    .Data = bb_fit,
    class = "pensynthboot"
  ))
}

#' @method predict pensynthboot
#' @export
predict.pensynthboot <- function(object, newdata, ...) {
  sapply(object, predict, newdata = newdata)
}


