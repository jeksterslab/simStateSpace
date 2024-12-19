#' Vectorize a Matrix
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a vector.
#'
#' @param x Matrix.
#'
#' @family Vectorization Functions
#' @keywords linearAlgebra vectorization internal
#' @noRd
.Vec <- function(x) {
  dim(x) <- NULL
  return(x)
}
