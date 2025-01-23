#' Half-Vectorize a Matrix
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Matrix.
#'
#' @return Returns a vector.
#'
#' @family Vectorization Functions
#' @keywords linearAlgebra vectorization internal
#' @noRd
.Vech <- function(x) {
  return(
    x[
      lower.tri(
        x = x,
        diag = TRUE
      )
    ]
  )
}
