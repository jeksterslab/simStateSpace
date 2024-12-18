#' Percentile Confidence Intervals
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param thetahatstar Numeric vector.
#'   Sampling distribution.
#' @param thetahat Numeric.
#'   Parameter estimate.
#' @param probs Numeric vector.
#'   Vector of probabilities corresponding to alpha level.
#'
#' @return Returns a matrix of estimates, standard errors,
#'   number of replications, and confidence intervals.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.PCCI <- function(thetahatstar,
                  thetahat,
                  probs) {
  thetahatstar <- as.vector(thetahatstar)
  thetahatstar <- thetahatstar[stats::complete.cases(thetahatstar)]
  ci <- stats::quantile(
    x = thetahatstar,
    probs = probs,
    names = FALSE
  )
  return(
    .CIFormat(
      thetahatstar = thetahatstar,
      thetahat = thetahat,
      probs = probs,
      ci = ci
    )
  )
}
