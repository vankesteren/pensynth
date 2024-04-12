#' Bayesian bootstrap inference
#'
#' Produce Bayesian bootstrap prediction intervals and ATE intervals
#' for pensynth objects. This is done through reweighing the columns
#' of the data (X0 and X1) using Dirichlet weights.
#'
#' @param object a fitted `pensynth` or `cvpensynth` object
#' @param n_iter `int` number of bootstrap iterations (default 1000)
#'
#' @export
bootstrap <- function(object, n_iter = 1000L) {
  UseMethod("bootstrap")
}

#' @rdname bootstrap
#' @export
bootstrap.pensynth <- function(object, n_iter = 1000L) {
  # create reweighed predictions
  bb_fit <- lapply(
    X = 1:n_iter,
    FUN = function(n) {
      # rerun original call but with boot weights
      cat("\rProgress: [", n, "/", n_iter, "]")
      update(object = object, boot_weight = TRUE)
    }
  )
  cat("\n")
  return(structure(
    .Data = list(boot_fit = bb_fit),
    class = "pensynthboot"
  ))
}

#' @rdname bootstrap
#' @export
bootstrap.cvpensynth <- function(object, n_iter = 1000L) {
  bootstrap.pensynth(object, n_iter)
}

#' @method predict pensynthboot
#' @export
predict.pensynthboot <- function(object, newdata, ...) {
  sapply(object$boot_fit, predict, newdata = newdata)
}


