#' Percentile Probabilities
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param alpha Numeric vector.
#'   Significance level.
#' @return Returns a vector of probabilities.
#'
#' @family Confidence Intervals Functions
#' @keywords nBootstrap ci internal
#' @noRd
.PCProbs <- function(alpha) {
  alpha <- sort(alpha)
  prob_ll <- alpha / 2
  prob_ul <- rev(1 - prob_ll)
  return(
    c(
      prob_ll,
      prob_ul
    )
  )
}
