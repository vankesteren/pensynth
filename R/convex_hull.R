#' Check whether treated unit is in the convex hull of donors
#'
#' This function finds out if the treated unit is in the convex hull
#' of the donor units.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param ... additional arguments passed to [clarabel::clarabel()]
#'
#' @details
#' This function does not actually construct the convex hull (which
#' is infeasible in higher dimensions), but rather it checks whether
#' the following linear program has a feasible solution:
#'
#' min q'w s.t. Aw = b
#'
#' with w constrained to be above 0 and sum to 1, the feasibility
#' of this linear program directly corresponds to whether the treated
#' is in the convex hull
#'
#' @return `bool` whether the treated unit is in the convex hull of
#' the donor units.
#'
#' @export
in_convex_hull <- function(X1, X0, ...) {
  # using linear progam feasibility
  # to find whether X1 is in
  # convex hull of X0.
  # see https://stackoverflow.com/a/43564754
  # clarabel solves q'w s.t. Aw = b

  # find dimensions of problem
  N <- ncol(X0)
  P <- nrow(X0)

  # First create the A matrix, but sparse for efficiency
  # A <- rbind(X0, 1, diag(N))
  A <- Matrix::sparseMatrix(
    i = c(rep(1:P, each = N), rep(P + 1, N), (P + 2):(P + 1 + N)),
    j = c(rep(1:N, P), 1:N, 1:N),
    x = c(t(X0), rep(1, N), rep(-1, N)),
    repr = "C"
  )

  # Then, create b
  b <- c(
    X1,       # X0w should equal X1
    1,        # sum(w) should equal 1
    rep(0, N) # -w[n] should be less than 0 for all N
  )

  # then create obj fun
  q <- rep(0, N) # we don't actually care about optim!

  # create cone list
  cone_list <- list(
    z = P + 1L, # we have P + 1 equality constraints
    l = N       # we have N inequality constraints
  )

  # now, just use clarabel to check for feasibility of solution
  sol <- clarabel::clarabel(A, b, q, cones = cone_list, ...)

  # if it was solved, return TRUE
  if (sol$status %in% c(2, 5)) return(TRUE)

  # if proven infeasible, return FALSE
  if (sol$status %in% c(3, 4, 6, 7)) return(FALSE)

  # else, return error with message
  stop(
    "Could not determine whether X0 is in convex hull:\n",
    clarabel::solver_status_descriptions()[sol$status]
  )
}
