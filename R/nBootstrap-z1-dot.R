#' z1
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param probs Numeric vector.
#'   Vector of probabilities corresponding to alpha level.
#'
#' @return Returns a numeric vector.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.Z1 <- function(probs) {
  return(
    stats::qnorm(p = probs)
  )
}
