#' Placebo permutation test for pensynth
#'
#' Perform a permutation test on a pensynth object, in the sense
#' of Abadie, Diamond, and Hainmueller (2010). The pensynth
#' method is performed multiple times, treating each donor as the
#' treated unit and the treated unit with the remaining donors as
#' the donor units.
#'
#' @param object a fitted `pensynth` or `cvpensynth` object
#' @param Y1 the post-intervention outcome of the treated unit
#' @param Y0 the post-intervention outcome of the donor units
#' (with N_donors columns)
#'
#'
#' @return A list with two elements
#' - E1, the treated unit effect, computed as `Y1 - Y0 %*% w`
#' - E0, the donor unit effects, computed in the same way but
#' using the permutation test's weights.
#' - ATE1, the estimated ATE of the treated unit
#' - ATE0, the estimated ATE of the donor units
#'
#' @references Abadie, A., Diamond, A., & Hainmueller, J. (2010).
#' Synthetic control methods for comparative case studies:
#' Estimating the effect of California’s tobacco control program.
#' Journal of the American statistical Association, 105(490),
#' 493-505.
#'
#' @export
placebo_test <- function(x, ...){
  UseMethod("placebo_test")
}

#' @rdname placebo_test
#' @export
placebo_test.pensynth <- function(object, Y1, Y0) {
  # original obs - pred
  est <- Y1 - predict(object, Y0)

  # create a list of weights
  N_donors <- length(object$w)
  ps_list <- lapply(
    1:N_donors,
    function(n) {
      # get X1 and X0 from original environment
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

  # compute test statistics
  ATE0 <- colMeans(null)
  ATE1 <- mean(est)

  return(structure(
    .Data = list(
      E0 = null,
      E1 = est,
      ATE1 = ATE1,
      ATE0 = ATE0
    ),
    class = "pensynthtest"
  ))
}

#' @rdname placebo_test
#' @export
placebo_test.cvpensynth <- function(object, Y1, Y0) {
  # original obs - pred
  est <- Y1 - predict(object, Y0)

  # create a list of weights
  N_donors <- length(object$w_opt)
  ps_list <- lapply(
    1:N_donors,
    function(n) {
      # get X1 and X0 from original environment
      new_X1 <- eval(getCall(object)$X0)[, n, drop = FALSE]
      new_Z1 <- eval(getCall(object)$Z0)[, n, drop = FALSE]
      new_X0 <- cbind(
        eval(getCall(object)$X1),
        eval(getCall(object)$X0)[, -n, drop = FALSE]
      )
      new_Z0 <- cbind(
        eval(getCall(object)$Z1),
        eval(getCall(object)$Z0)[, -n, drop = FALSE]
      )
      update(object = object, X1 = new_X1, X0 = new_X0, Z1 = new_Z1, Z0 = new_Z0)
    }
  )

  # create prediction for each iter
  null <- sapply(
    1:N_donors,
    function(n) {
      Y0[, n, drop = FALSE] - predict(ps_list[[n]], cbind(Y1, Y0[, -n, drop = FALSE]))
    }
  )

  # compute test statistics
  ATE0 <- colMeans(null)
  ATE1 <- mean(est)

  return(structure(
    .Data = list(
      E0 = null,
      E1 = est,
      ATE1 = ATE1,
      ATE0 = ATE0
    ),
    class = "pensynthtest"
  ))
}

#' Plotting a pensynth permutation object
#'
#' Plotting the reference distribution and the
#' estimated treatement effect for the treated unit
#' for the pensynth permutation test.
#'
#' @param object a `pensynthtest` object
#' @param ... additional parameters passed to `plot`
#'
#' @seealso [base::plot()]
#'
#' @method plot pensynthtest
#'
#' @export
plot.pensynthtest <- function(object, ...) {
  val_range <- range(c(object$E0, object$E1))
  N_post <- nrow(object$E0)
  N_donors <- ncol(object$E0)
  plot(
    NA, ylim = val_range, xlim = c(1, N_post),
    ylab = "Treatment effect",
    xlab = "Post-intervention timepoint",
    ...
  )
  for (n in 1:N_donors) lines(object$E0[,n], col = "grey")
  lines(object$E1, lwd = 1.2)
  legend("topright", lty = c(1, 1), col = c("grey", "black"),
         legend = c("reference", "estimate"))
}
