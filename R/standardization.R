#' Function to scale the X matrices in synthetic control
#'
#' The default implementation of synthetic control first
#' scales the X matrices using the mean and standard
#' deviation of the X1 and X0 matrices together.
#'
#' @param X1 the X1 matrix
#' @param X0 the X0 matrix
#'
#' @details See the pre-processing in [Synth::synth()].
#'
#' @seealso [base::scale()]
#'
#' @returns a list with X1 and X0 standardized.
#'
#' @importFrom stats sd
#'
#' @keywords internal
standardize_X <- function(X1, X0) {
  X <- cbind(X1, X0)
  mu <- apply(X, 1, mean)
  sc <- apply(X, 1, sd)
  return(list(X1 = (X1-mu)/sc, X0 = (X0-mu)/sc))
}
