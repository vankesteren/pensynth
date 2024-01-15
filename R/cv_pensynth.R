#' Cross-validated penalized synthetic control estimator
#'
#' Compute a penalized synthetic control estimator with cross-validation for the
#' lambda penalty parameter. Lambda will be determined by minimizing the mean squared
#' error on a hold-out set of pre-intervention outcome time-series.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param v `N_covars vector` of variable weights
#' @param Z1 `N_targets by 1 matrix` of treated unit hold-out outcome
#' @param Z0 `N_targets by N_donors matrix` of donor unit hold-out outcome
#' @param nlambda `integer` length of lambda sequence (see details)
#' @param opt_pars `osqp` settings using [osqp::osqpSettings()]
#'
#' @details The lambda sequence is an exponentially increasing sequence where
#' The minimum lambda is always 1e-7, the max lambda is determined by the data.
#'
#' @return A list of the lambda sequence, the associated weights, and the mses
#'
#' @seealso [pensynth()] [plot_path()]
#'
#' @importFrom utils capture.output
#'
#' @examples
#' set.seed(45)
#' N_covar <- 7
#' N_donor <- 50
#' N_target <- 12
#'
#' w  <- runif(N_donor)
#' w[5:N_donor] <- 0
#' w  <- w / sum(w)
#' v  <- rep(1, N_covar)
#' X0 <- matrix(rnorm(N_covar*N_donor), N_covar)
#' X1 <- X0%*%w
#' Z0 <- matrix(rnorm(N_target*N_donor), N_target)
#' Z1 <- Z0%*%w
#'
#' res <- cv_pensynth(X1, X0, v, Z1, Z0)
#' plot_path(res)
#'
#' @export
cv_pensynth <- function(X1, X0, v, Z1, Z0, nlambda = 100, opt_pars = osqp::osqpSettings(polish = TRUE)) {
  N_donors <- ncol(X0)
  X0v <- X0*sqrt(v)
  X1v <- X1*sqrt(v)

  X0VX0 <- crossprod(X0v)
  X1VX0 <- crossprod(X1v, X0v)
  Delta <- apply(X0v - c(X1v), 2, crossprod)

  lseq <- lambda_sequence(X1VX0, Delta, nlambda)

  # linear constraint matrix
  Amat   <- rbind(rep(1, N_donors), diag(N_donors))
  lbound <- c(1, rep(0, N_donors))
  ubound <- rep(1, N_donors + 1)

  w_path <- sapply(lseq, function(lambda) {
    o <- capture.output({
      # Instantiate the quadratic program solver
      qpsolver <- osqp::osqp(
        P = X0VX0,
        q = -X1VX0 + lambda*Delta,
        A = Amat,
        l = lbound,
        u = ubound,
        pars = opt_pars
      )
      # solve
      result <- qpsolver$Solve()
    })
    return(result$x)
  })
  colnames(w_path) <- lseq
  e_path <- sapply(1:nlambda, \(i) crossprod(Z1 - Z0 %*% w_path[,i])) / length(Z1)

  out_obj <- structure(
    .Data = list(
      w_opt    = w_path[,which.min(e_path)],
      l_opt    = lseq[which.min(e_path)],
      lseq     = lseq,
      w_path   = w_path,
      mse_path = e_path
    ),
    class = "cvpensynth"
  )

  return(out_obj)
}

#' Plotting for cross-validated penalized synthetic control objects
#'
#' Displays a mean squared error curve and weights curve as a function
#' of lambda, the penalization parameter.
#'
#' @param object a `cvpensynth` output object
#' @param ... additional arguments passed to `plot()`
#'
#' @seealso [cv_pensynth()] [pensynth()]
#'
#' @importFrom graphics lines par
#'
#' @export
plot_path <- function(object, ...) {
  nw <- nrow(object$w_path)
  mfrow_old <- par("mfrow")
  on.exit(par(mfrow = mfrow_old))
  par(mfrow = c(2, 1))
  plot(
    object$lseq,
    object$mse_path,
    log = "x",
    ylab = "MSE",
    xlab = "Lambda",
    type = "l",
    main = "Mean squared prediction errors",
    ...
  )
  plot(
    object$lseq,
    object$w_path[1, ],
    log = "x",
    ylab = "Weight",
    xlab = "Lambda",
    type = "l",
    ylim = c(0, 1),
    main = "Weights",
    ...
  )
  for (i in 2:nw) {
    lines(object$lseq, object$w_path[i, ], lty = i)
  }
}

lambda_sequence <- function(X1VX0, Delta, nlambda) {
  lmin <- 1e-7
  lmax <- sum(abs(X1VX0/Delta))
  return(exp(seq(log(lmin), log(lmax), len = nlambda)))
}
