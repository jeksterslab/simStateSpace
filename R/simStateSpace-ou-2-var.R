#' Convert Parameters from the Ornstein–Uhlenbeck Model
#' to the Vector Autoregressive Model
#'
#' This function converts parameters from the Ornstein–Uhlenbeck model
#' to the vector autoregressive model parameterization.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param delta_t Positive integer.
#'   Discrete time interval (\eqn{\delta_t}).
#' @inheritParams SimSSMOU
#'
#' @return Returns a list of vector autoregressive model parameters:
#'   - `alpha`: Numeric vector.
#'     Vector of intercepts for the dynamic model
#'     (\eqn{\boldsymbol{\alpha}}).
#'   - `beta`: Numeric matrix.
#'     Transition matrix relating the values of the latent variables
#'     at time `t - delta_t` to those at time `t`
#'     where `delta_t` is a positive integer.
#'     (\eqn{\boldsymbol{\beta}}).
#'   - `psi`: Numeric matrix.
#'     The process noise covariance matrix
#'     (\eqn{\boldsymbol{\Psi}}).
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
#' p <- k <- 3
#' mu <- c(0, 0, 0)
#' phi <- matrix(
#'   data = c(-0.357, 0.771, -0.450, 0, -0.511, 0.729, 0, 0, -0.693),
#'   nrow = p
#' )
#' sigma <- diag(p)
#' delta_t <- 1
#'
#' OU2VAR(
#'   mu = mu,
#'   phi = phi,
#'   sigma = sigma,
#'   delta_t = delta_t
#' )
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
