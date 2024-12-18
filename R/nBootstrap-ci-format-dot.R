#' Format Confidence Intervals
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param thetahatstar Numeric vector.
#'   Sampling distribution.
#' @param thetahat Numeric.
#'   Parameter estimate.
#' @param probs Numeric vector.
#'   Vector of probabilities corresponding to alpha level.
#' @param ci Numeric vector.
#'   Confidence intervals.
#'
#' @return Returns a vector of
#'   estimate,
#'   standard error of estimate,
#'   number of replications,
#'   and
#'   confidence intervals.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.CIFormat <- function(thetahatstar,
                      thetahat,
                      probs,
                      ci) {
  out <- c(
    thetahat,
    stats::sd(thetahatstar),
    length(thetahatstar),
    ci
  )
  names(out) <- c(
    "est",
    "se",
    "R",
    paste0(
      probs * 100,
      "%"
    )
  )
  return(out)
}
