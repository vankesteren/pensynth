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
#' @param return_solver_info `boolean` whether to return diagnostic information concerning solver (default FALSE)
#' @param verbose whether to print progress messages. Default on if in an interactive session.
#'
#' @details The lambda sequence is an exponentially increasing sequence where
#' The minimum lambda is always 1e-11, the max lambda is determined by the data.
#'
#' @returns A list of the lambda sequence, the associated weights, and the mses. If
#' `return_solver_info` is `TRUE`, the list will also contain diagnostic information about
#' the solvers.
#'
#' @seealso [pensynth()], [plot.cvpensynth()], [placebo_test()], [simulate_data()]
#'
#' @importFrom utils capture.output
#' @importFrom stats loess
#'
#' @example R/examples/example_cv_pensynth.R
#'
#' @export
cv_pensynth <- function(
  X1, X0, Z1, Z0, v = 1, nlambda = 100,
  opt_pars = clarabel::clarabel_control(),
  standardize = TRUE, return_solver_info = FALSE,
  verbose = interactive(), adaptive_lambda = TRUE
) {

  N_treated <- if (is.matrix(X1)) ncol(X1) else 1
  if (N_treated > 1) {
    if (verbose) cli::cli_progress_step("{N_treated} treated units detected, estimating multiple models")

    # a little recursion :)
    results <- vector("list", N_treated)
    if (verbose) pb <- cli::cli_progress_bar("Fitting models.", total = N_treated, current = FALSE)
    for (n in 1:N_treated) {
      results[[n]] <- cv_pensynth(
        X1 = X1[, n, drop=FALSE],
        X0 = X0,
        Z1 = Z1[, n, drop=FALSE],
        Z0 = Z0,
        v = v,
        nlambda = nlambda,
        opt_pars = opt_pars,
        standardize = standardize,
        return_solver_info = return_solver_info,
        verbose = FALSE
      )
      if (verbose) cli::cli_progress_update(id = pb, force = TRUE)
    }

    if (!adaptive_lambda) {
      # compute a shared lambda
      if (verbose) cli::cli_progress_step("Estimating common penalty parameter.")
      lams <- c(sapply(results, function(x) x$lseq))
      mses <- c(sapply(results, function(x) x$mse_path))

      # Flexible nonparametric local regression model to estimate the minimum
      mse_mod <- loess(sqrt(mses) ~ log(lams), span = 0.05) # loess with 5% neighbourhood
      mse_est <- predict(mse_mod)^2
      lam_opt <- lams[which.min(mse_est)]

      # Extract weights matrix
      W_opt <- sapply(results, function(x) {
        opt_idx <- which.min(abs(x$lseq - lam_opt))
        x$w_path[, opt_idx]
      })
    } else {
      # Use individual lambda per treated unit
      W_opt   <- sapply(results, function(x) x$w_opt)
      lam_opt <- vapply(results, function(x) x$l_opt, 0.0)
    }

    if (verbose) cli::cli_progress_step("Collecting output.")

    # Construct a list of outputs
    out_obj <- list(
      w_opt    = W_opt,
      l_opt    = lam_opt,
      lseq     = lapply(results, \(x) x$lseq),
      w_path   = lapply(results, \(x) x$w_path),
      mse_path = lapply(results, \(x) x$mse_path),
      call     = match.call()
    )

    # Convert the list to a cvpensynth object
    out_obj <- structure(
      .Data = out_obj,
      class = "cvpensynth"
    )
    return(out_obj)
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
  solve_qp <- function(id) {
    lambda <- lseq[id]
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

  if (verbose) {
    cli::cli_progress_step("Fitting models.", name = "fitmsg")
    solver_output <- sapply(
      X = cli::cli_progress_along(x = lseq, name = "  Fitting models.", current = FALSE),
      FUN = solve_qp
    )
  } else {
    solver_output <- sapply(X = seq_along(lseq), FUN = solve_qp)
  }

  if (verbose) cli::cli_progress_step("Getting output.")

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
  nt <- if (is.matrix(x$w_opt)) ncol(x$w_opt) else 1
  nd <- if (is.matrix(x$w_opt)) nrow(x$w_opt) else length(x$w_opt)
  w8 <- if (nt > 1) x$w_opt[1:min(nd, 8), 1] else x$w_opt[1:min(nd, 8)]
  cat("Hold-out validated pensynth model\n---------------------------------\n")
  print(x$call)
  cat("\n- Treated units:", nt, "\n")
  cat("- Donor units:", nd, "\n")
  cat("- lambda:", x$l_opt, "\n")
  cat("- mse:", get_mse_cvpensynth(x), "\n")
  cat("- w:", round(w8, 3))
  if(length(x$w_opt) > 8) cat("...")

  return(invisible(x))
}

#' Get average mse for cvpensynth object
#'
#' @param fit The cvpensynth object
#'
#' @returns numeric
#'
#' @keywords internal
get_mse_cvpensynth <- function(fit) {
  # sorry for this horrible nested list code
  # this computes average MSE over the treated unit
  # models
  if (ncol(x$w_opt) == 1) return(round(min(x$mse_path), 3))

  if (isFALSE(x$call[["adaptive_lambda"]])) return(
    mean(sapply(seq_along(x$mse_path), function(i) {
      opt_idx <- which.min(abs(x$lseq[[i]] - x$l_opt))
      x$mse_path[[i]][opt_idx]
    }))
  )
  return(round(mean(sapply(x$mse_path, min)), 3))
}


#' Plotting for hold-out validated penalized synthetic control objects
#'
#' Displays a mean squared error curve and weights curve as a function
#' of lambda, the penalization parameter.
#'
#' @param x a `cvpensynth` output object
#' @param treated_unit index of the treated unit to display
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
plot.cvpensynth <- function(x, treated_unit = 1, ...) {

  # collect info based on number of treated units
  n_trt <- ncol(x$w_opt)
  if (n_trt == 1) {
    lseq     <- x$lseq
    mse_path <- x$mse_path
    w_path   <- x$w_path
    l_opt    <- x$l_opt
  } else {
    single_lambda <- isFALSE(x$call[["adaptive_lambda"]])

    lseq     <- x$lseq[[treated_unit]]
    mse_path <- x$mse_path[[treated_unit]]
    w_path   <- x$w_path[[treated_unit]]
    l_opt    <- if (single_lambda) x$l_opt else x$l_opt[treated_unit]
  }

  # start plotting
  mfrow_old <- par("mfrow")
  on.exit(par(mfrow = mfrow_old))
  par(mfrow = c(2, 1))

  nw <- nrow(w_path)
  plot(
    x    = lseq,
    y    = mse_path,
    log  = "x",
    ylab = "MSE",
    xlab = "Lambda",
    type = "l",
    main = "Mean squared prediction errors",
    ...
  )
  abline(v = l_opt, col = "grey")
  plot(
    x    = lseq,
    y    = w_path[1, ],
    log  = "x",
    ylab = "Weight",
    xlab = "Lambda",
    type = "l",
    ylim = c(0, 1),
    main = "Unit weights",
    ...
  )
  abline(v = l_opt, col = "grey")
  for (i in 2:nw) {
    lines(x = lseq, y = w_path[i, ], lty = i)
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
