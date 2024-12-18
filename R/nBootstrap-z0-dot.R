#' z0
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param thetahatstar Numeric vector.
#'   Sampling distribution.
#' @param thetahat Numeric.
#'   Parameter estimate.
#'
#' @return Returns a numeric vector of length one.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.Z0 <- function(thetahatstar,
                thetahat) {
  return(
    stats::qnorm(
      p = sum(
        thetahatstar < thetahat
      ) / length(thetahatstar)
    )
  )
}
