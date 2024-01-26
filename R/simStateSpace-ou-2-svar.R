#' Convert Parameters from the Ornstein–Uhlenbeck Model
#' to the Structural Vector Autoregressive Model
#'
#' This function converts parameters from the Ornstein–Uhlenbeck model
#' to the structural vector autoregressive model parameterization.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams OU2VAR
#'
#' @return Returns a list of structural vector autoregressive model parameters:
#'   - `alphastar`: Numeric vector.
#'     Vector of intercepts for the dynamic model
#'     (\eqn{\boldsymbol{\alpha}^{\ast}}).
#'   - `betastar`: Numeric matrix.
#'     Transition matrix relating the values of the latent variables
#'     at time `t - delta_t` to those at time `t`
#'     where `delta_t` is a positive integer.
#'     (\eqn{\boldsymbol{\beta}^{\ast}}).
#'   - `psistar`: Numeric matrix.
#'     The process noise covariance matrix
#'     (\eqn{\boldsymbol{\Psi}^{\ast}}).
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
#' OU2SVAR(
#'   mu = mu,
#'   phi = phi,
#'   sigma = sigma,
#'   delta_t = delta_t
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ou
#' @export
OU2SVAR <- function(mu,
                    phi,
                    sigma,
                    delta_t) {
  delta_t <- as.integer(delta_t)
  var <- OU2SSM(
    mu = mu,
    phi = phi,
    sigma = sigma,
    delta_t = delta_t
  )
  p <- length(mu)
  iden <- diag(p)
  alphastar <- -1 * (
    (
      (-1)^(p + 1)
    ) * (
      delta_t^(p)
    ) * (
      solve(var$beta)
    ) - iden
  )
  bread <- iden - alphastar
  out <- list(
    alphastar = alphastar,
    betastar = bread %*% var$beta,
    psistar = bread %*% var$psi %*% t(bread)
  )
  return(out)
}
