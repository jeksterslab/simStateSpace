#' Steady-State Covariance Matrix for the
#' Latent Variables in the
#' State Space Model
#'
#' The steady-state covariance matrix
#' for the latent variables
#' in the state space model
#' \eqn{\mathrm{Cov} \left( \boldsymbol{\eta} \right)}
#' is given by
#' \deqn{
#'   \mathrm{vec}
#'   \left(
#'     \mathrm{Cov} \left( \boldsymbol{\eta} \right)
#'   \right)
#'   =
#'   \left(
#'     \mathbf{I} - \boldsymbol{\beta} \otimes \boldsymbol{\beta}
#'   \right)^{-1}
#'   \mathrm{vec} \left( \boldsymbol{\Psi} \right)
#' }
#' where
#' \eqn{\boldsymbol{\beta}}
#' is the transition matrix relating the values of the latent variables
#' at the previous to the current time point
#' and
#' \eqn{\boldsymbol{\Psi}}
#' is the covariance matrix
#' of volatility or randomness in the process.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param psi Numeric matrix.
#'   The covariance matrix
#'   of the process noise
#'   (\eqn{\boldsymbol{\Psi}}).
#' @inheritParams SimSSMFixed
#'
#' @examples
#' beta <- 0.50 * diag(3)
#' psi <- 0.001 * diag(3)
#' SSMCovEta(
#'   beta = beta,
#'   psi = psi
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SSMCovEta <- function(beta, psi) {
  .SSMCovEta(
    beta = beta,
    psi = psi
  )
}
