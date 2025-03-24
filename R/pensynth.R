#' Penalized synthetic control estimator
#'
#' For a given set of variable weights (v) this function estimates
#' the unit weights for a synthetic control with penalization
#' according to Abadie & L'Hour (2021). This function deals with only a
#' single treated unit.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param v `N_covars vector` of variable weights (default 1)
#' @param lambda `numeric` penalization parameter
#' @param opt_pars `clarabel` settings using [clarabel::clarabel_control()]
#' @param standardize `boolean` whether to standardize the input matrices (default TRUE)
#' @param verbose whether to print progress messages. Default on if in an interactive session.
#'
#' @details This routine uses the same notation of the original [Synth::synth()] implementation
#' but uses a different, faster quadratic program solver (namely, [clarabel::clarabel()]).
#' Additionally, it implements the penalization procedure described in Abadie & L'Hour (2021),
#' such that the loss function is as in equation 5 of that paper (but for a single treated unit).
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
#' @returns A list with two values: `w`, the estimated weights; and
#' `solution`, the result of the optimization.
#'
#' @importFrom utils capture.output
#'
#' @example R/examples/example_pensynth.R
#'
#' @seealso [cv_pensynth()], [placebo_test()], [simulate_data()], [Synth::synth()]
#'
#' @export
pensynth <- function(X1,
                     X0,
                     v = 1,
                     lambda = 0,
                     opt_pars = clarabel::clarabel_control(),
                     standardize = TRUE,
                     verbose = interactive()) {

  N_treated <- ncol(X1)
  if (N_treated > 1) {
    if (verbose) cli::cli_alert_info("Multiple treated units detected, estimating multiple models")

    # a little recursion :)
    results <- vector("list", N_treated)
    if (verbose) cli::cli_progress_bar("Fitting models.")
    for (n in 1:N_treated) {
      results[[n]] <- pensynth(X1[,n,drop=FALSE], X0, v, lambda, opt_pars, standardize, verbose = FALSE)
      if (verbose) cli::cli_progress_update()
    }

    return(structure(
      .Data = list(
        w = sapply(results, \(x) x[["w"]]),
        solution = lapply(results, \(x) x[["result"]]),
        call = match.call()
      ),
      class = "pensynth"
    ))
  }

  if (verbose) cli::cli_progress_step("Preparing data.")
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

  # Constraint matrices (sparse for efficiency)
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

  # Run the quadratic program solver
  if (verbose) cli::cli_progress_step("Fitting model.")
  result <- clarabel::clarabel(
    P = X0VX0,
    q = -X1VX0 + lambda*Delta,
    A = Amat,
    b = B,
    cones = list(
      z = 1L, # There is 1 equality
      l = N_donors # There are N_donors inequalities
    ),
    control = opt_pars,
  )

  # clarabel only returns a numeric status code, so we'll add a
  # human-readable status column here (plus a description)
  result$status_description <- clarabel::solver_status_descriptions()[result$status][[1]]
  result$status <- names(clarabel::solver_status_descriptions()[result$status])


  return(structure(
    .Data = list(
      w = result$x,
      solution = result,
      call = match.call()
    ),
    class = "pensynth"
  ))
}

#' Print pensynth model
#'
#' @param x a pensynth object
#' @param ... ignored
#'
#' @method print pensynth
#'
#' @returns the pensynth object, invisibly
#'
#' @export
print.pensynth <- function(x, ...) {
  cat("Pensynth model\n--------------\n")
  cat("- call: ")
  print(x$call)
  cat("- solution:", x$solution$status, "\n")
  cat("- w:", round(x$w, 3)[1:min(length(x$w), 8)])
  if(length(x$w) > 8) cat("...")
  return(invisible(x))
}


#' Create prediction from pensynth model
#'
#' Matrix multiplies the values in `newdata` by the unit weights
#' extracted from the pensynth object to produce predicted
#' values.
#'
#' @param object a fitted cvpensynth model
#' @param newdata N_values * N_donors matrix of
#' values for the donor units.
#' @param ... ignored
#'
#' @returns a matrix (column vector) of predicted values
#'
#' @importFrom stats predict approx
#'
#' @method predict pensynth
#'
#' @export
predict.pensynth <- function(object, newdata, ...) {
  return(newdata %*% object$w)
}
