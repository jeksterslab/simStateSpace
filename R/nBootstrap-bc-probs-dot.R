#' Bias Corrected Probablities
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param z0 Numeric.
#' @param z1 Numeric vector
#'   with length equal to two times the length of alpha.
#' @return Returns a vector of prababilities.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.BCProbs <- function(z0,
                     z1) {
  return(
    stats::pnorm(
      q = 2 * z0 + z1
    )
  )
}
