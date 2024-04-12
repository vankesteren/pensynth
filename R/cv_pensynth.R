#' Hold-out validated penalized synthetic control estimator
#'
#' Compute a penalized synthetic control estimator with hold-out validation for the
#' lambda penalty parameter. Lambda will be determined by minimizing the mean squared
#' error on a hold-out set of pre-intervention outcome time-series.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param Z1 `N_targets by 1 matrix` of treated unit hold-out outcome
#' @param Z0 `N_targets by N_donors matrix` of donor unit hold-out outcome
#' @param v `N_covars vector` of variable weights, default 1
#' @param nlambda `integer` length of lambda sequence (see details)
#' @param opt_pars `clarabel` settings using [clarabel::clarabel_control()]
#' @param standardize `boolean` whether to standardize the input matrices (default TRUE)
#' @param boot_weight `boolean` whether to apply Bayesian bootstrap weights (see details)
#' @param return_solver_info `boolean` whether to return diagnostic information concerning solver (default FALSE)
#'
#' @details The lambda sequence is an exponentially increasing sequence where
#' The minimum lambda is always 1e-11, the max lambda is determined by the data.
#'
#' Bayesian bootstrap weights can be applied as well, to produce uncertainty estimates for the
#' predictions. The weights will be generated randomly on-the-fly.
#'
#' @returns A list of the lambda sequence, the associated weights, and the mses. If
#' `return_solver_info` is `TRUE`, the list will also contain diagnostic information about
#' the solvers.
#'
#' @seealso [pensynth()], [plot.cvpensynth()], [placebo_test()], [simulate_data()]
#'
#' @importFrom utils capture.output
#'
#' @example R/examples/example_cv_pensynth.R
#'
#' @export
cv_pensynth <- function(X1, X0, Z1, Z0, v = 1, nlambda = 100, opt_pars = clarabel::clarabel_control(),
                        standardize = TRUE, boot_weight = FALSE, return_solver_info = FALSE) {
  if (standardize) {
    st <- standardize_X(X1, X0)
    X0 <- st$X0
    X1 <- st$X1
  }
  N_donors <- ncol(X0)
  X0v <- X0*sqrt(v)
  X1v <- X1*sqrt(v)

  if (boot_weight) {
    N_covar <- nrow(X0)
    N_pre <- nrow(Z0)
    w <- rexp(N_donors)
    X0v <- X0v * rep(sqrt(w), each = N_covar)
    Z0  <- Z0 * rep(sqrt(w), each = N_pre)
    sw1 <- sqrt(rexp(1))
    X1v <- X1v * sw1
    Z1  <- Z1 * sw1
  }

  X0VX0 <- crossprod(X0v)
  X1VX0 <- crossprod(X1v, X0v)
  Delta <- apply(X0v - c(X1v), 2, crossprod)

  lseq <- lambda_sequence(X1VX0, Delta, nlambda)

  # Constraint matrices
  # Amat <- rbind(
  #   rep(1, N_donors), # Sum to 1 constraint
  #   -diag(N_donors) # Individ. weights gte 0 constraint
  # )
  Amat <- Matrix::sparseMatrix(
    i = c(rep(1, N_donors), 2:(N_donors + 1)),
    j = c(1:N_donors, 1:N_donors),
    x = c(rep(1, N_donors), rep(-1, N_donors)),
    repr = "C"
  )
  B <- c(
    1, # Sum to 1 constraint
    rep(0, N_donors) # Individ. weights gte 0 constraint
  )

  # Define function for solving qp for a given lambda
  solve_qp <- function(lambda) {
    # run the quadratic program solver
    result <- clarabel::clarabel(
      P = X0VX0,
      q = -X1VX0 + lambda*Delta,
      A = Amat,
      b = B,
      cones = list(
        z = 1L, # There is 1 equality
        l = N_donors # There are N_donors inequalities
      ),
      control = opt_pars
    )

    # clarabel only returns a numeric status code, so we'll add a
    # human-readable status column here (plus a description)
    result$status_description <- clarabel::solver_status_descriptions()[result$status][[1]]
    result$status <- names(clarabel::solver_status_descriptions()[result$status])

    # Return result
    return(result)
  }

  solver_output <- sapply(lseq, solve_qp)

  # Extract weights
  w_path <- do.call(cbind, solver_output["x", ])

  colnames(w_path) <- lseq
  e_path <- sapply(1:nlambda, \(i) crossprod(Z1 - Z0 %*% w_path[,i])) / length(Z1)

  # Construct a list of outputs
  out_obj <- list(
      w_opt    = w_path[,which.min(e_path)],
      l_opt    = lseq[which.min(e_path)],
      lseq     = lseq,
      w_path   = w_path,
      mse_path = e_path,
      call     = match.call()
  )

  # If we've been requested to return info about the solving process, do so
  if (return_solver_info) {
    # Remove unneeded columns from the solver output matrix
    rows_to_drop <- c("x", "y", "s", "z")
    solver_output <- solver_output[!rownames(solver_output) %in% rows_to_drop, ]

    # Add each row from the solver output matrix to .Data
    for (i in 1:nrow(solver_output)) {
      row_name <- rownames(solver_output)[i]
      out_obj[[row_name]] <- unlist(solver_output[i, ])
    }
  }

  # Convert the list to a cvpensynth object
  out_obj <- structure(
    .Data = out_obj,
    class = "cvpensynth"
  )

  return(out_obj)
}

#' Print cvpensynth model
#'
#' @param x a cvpensynth object
#' @param ... ignored
#'
#' @method print cvpensynth
#'
#' @returns the cvpensynth object, invisibly
#'
#' @export
print.cvpensynth <- function(x, ...) {
  cat("Hold-out validated pensynth model\n---------------------------------\n")
  cat("- call: ")
  print(x$call)
  cat("- lambda:", x$l_opt, "\n")
  cat("- mse:", round(min(x$mse), 3), "\n")
  cat("- w:", round(x$w_opt, 3)[1:min(length(x$w_opt), 8)])
  if(length(x$w_opt) > 8) cat("...")

  return(invisible(x))
}

#' Plotting for hold-out validated penalized synthetic control objects
#'
#' Displays a mean squared error curve and weights curve as a function
#' of lambda, the penalization parameter.
#'
#' @param x a `cvpensynth` output object
#' @param ... additional arguments passed to `plot()`
#'
#' @returns No return value, called for side effects
#'
#' @seealso [cv_pensynth()] [pensynth()]
#'
#' @importFrom graphics lines par abline
#'
#' @method plot cvpensynth
#'
#' @export
plot.cvpensynth <- function(x, ...) {
  nw <- nrow(x$w_path)
  mfrow_old <- par("mfrow")
  on.exit(par(mfrow = mfrow_old))
  par(mfrow = c(2, 1))
  plot(
    x$lseq,
    x$mse_path,
    log = "x",
    ylab = "MSE",
    xlab = "Lambda",
    type = "l",
    main = "Mean squared prediction errors",
    ...
  )
  abline(v = x$l_opt, col = "grey")
  plot(
    x$lseq,
    x$w_path[1, ],
    log = "x",
    ylab = "Weight",
    xlab = "Lambda",
    type = "l",
    ylim = c(0, 1),
    main = "Unit weights",
    ...
  )
  abline(v = x$l_opt, col = "grey")
  for (i in 2:nw) {
    lines(x$lseq, x$w_path[i, ], lty = i)
  }
}

#' Determine lambda sequence
#'
#' This function uses the weighted cross-product matrix
#' (X1VX0) and Delta to determine the lambda sequence.
#' This sequence will be exponentially increasing so it
#' is easy to plot with a logarithmic x-axis
#'
#' @param X1VX0 the weighted cross-product matrix
#' @param Delta the matching penalty matrix
#' @param nlambda the number of lambda values
#'
#' @details
#' The formula for the maximum lambda value was determined
#' empirically, with an eye for the form of the quadratic
#' program. In general, the max lambda should be so large
#' that we are practically in "nearest-neighbour" matching
#' territory. This formula ensures this for a wide range
#' of input parameters.
#'
#' @seealso [plot.cvpensynth()]
#'
#' @returns lambda sequence as a numeric vector
#'
#' @keywords internal
lambda_sequence <- function(X1VX0, Delta, nlambda) {
  lmin <- 1e-11
  lmax <- sum(abs(X1VX0/Delta))
  return(exp(seq(log(lmin), log(lmax), len = nlambda)))
}


#' Create prediction from cvpensynth model
#'
#' Matrix multiplies the values in `newdata` by the unit weights
#' extracted from the cvpensynth object to produce predicted
#' values.
#'
#' @param object a fitted cvpensynth model
#' @param newdata N_values * N_donors matrix of
#' values for the donor units.
#' @param lambda desired lambda value (defaults to optimal lambda)
#' @param ... ignored
#'
#' @details
#' For a chosen lambda that is not in the list of tested lambdas
#' in the cvpensynth object, the closest lambda (on the log scale)
#' will be chosen.
#'
#' @returns a matrix (column vector) of predicted values
#'
#' @importFrom stats predict approx
#'
#' @method predict cvpensynth
#'
#' @export
predict.cvpensynth <- function(object, newdata, lambda, ...) {
  if (missing(lambda)) return(newdata %*% object$w_opt)
  # find lambda idx
  lambda_idx <- which.min(abs(log(object[["lseq"]]) - log(lambda)))
  message("Closest lambda: ", object[["lseq"]][lambda_idx])
  return(newdata %*% object[["w_path"]][,lambda_idx])
}

#' @importFrom stats weights
#' @export
weights.cvpensynth <- function(object, ...) object$w_opt
