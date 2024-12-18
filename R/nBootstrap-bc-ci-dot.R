#' Bias Corrected Confidence Intervals
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param thetahatstar Numeric vector.
#'   Sampling distribution.
#' @param thetahat Numeric.
#'   Parameter estimate.
#' @param probs Numeric vector.
#'   Vector of probabilities.
#' @param z0 Numeric.
#' @param z1 Numeric vector
#'   with length equal to two times the length of alpha.
#' @return Returns a matrix of estimates, standard errors,
#'   number of replications, and confidence intervals.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.BCCI <- function(thetahatstar,
                  thetahat,
                  probs,
                  z0,
                  z1) {
  thetahatstar <- as.vector(thetahatstar)
  thetahatstar <- thetahatstar[stats::complete.cases(thetahatstar)]
  ci <- stats::quantile(
    x = thetahatstar,
    probs = .BCProbs(
      z0 = z0,
      z1 = z1
    ),
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
