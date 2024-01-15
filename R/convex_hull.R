#' Is the treated unit in the convex hull
#'
#' This function finds out if the treated unit is in the convex hull
#' of the donor units.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#'
#' @importFrom geometry convhulln inhulln
#'
#' @return `bool` whether the treated unit is in the convex hull of
#' the donor units.
#'
#' @seealso [geometry::convhulln()] [geometry::inhulln()]
#' @export
in_convex_hull <- function(X1, X0) {
  return(inhulln(convhulln(t(X0)), t(X1)))
}
