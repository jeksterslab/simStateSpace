#' Convert Parameters from the Ornstein–Uhlenbeck Model
#' to the Vector Autoregressive Model
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param delta_t Numeric.
#'   Discrete time interval (\eqn{\delta_t}).
#' @inheritParams OU2SSM
#'
#' @examples
#' p <- k <- 2
#' mu <- c(5.76, 5.18)
#' phi <- matrix(
#'   data = c(0.10, -0.05, -0.05, 0.10),
#'   nrow = p
#' )
#' sigma <- matrix(
#'   data = c(2.79, 0.06, 0.06, 3.27),
#'   nrow = p
#' )
#' delta_t <- 1
#'
#' OU2VAR(
#'   mu = mu,
#'   phi = phi,
#'   sigma = sigma,
#'   delta_t = delta_t
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ou
#' @export
OU2VAR <- function(mu,
                   phi,
                   sigma,
                   delta_t) {
  delta_t <- as.integer(delta_t)
  return(
    OU2SSM(
      mu = mu,
      phi = phi,
      sigma = sigma,
      delta_t = delta_t
    )
  )
}
