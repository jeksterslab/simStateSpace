#' Steady-State Covariance Matrix for the
#' Observed Variables in the
#' State Space Model
#'
#' The steady-state covariance matrix
#' for the observed variables
#' in the state space model
#' \eqn{\mathrm{Cov} \left( \mathbf{y} \right)}
#' is given by
#' \deqn{
#'   \mathrm{Cov} \left( \mathbf{y} \right)
#'   =
#'   \boldsymbol{\Lambda}
#'   \mathrm{Cov} \left( \boldsymbol{\eta} \right)
#'   \boldsymbol{\Lambda}^{\prime}
#'   +
#'   \boldsymbol{\Theta}
#' }
#' where
#' \eqn{\boldsymbol{\Lambda}}
#' is the matrix of factor loadings,
#' \eqn{\boldsymbol{\Theta}}
#' is the covariance matrix of
#' the measurement error,
#' and
#' \eqn{\mathrm{Cov} \left( \boldsymbol{\eta} \right)}
#' is the steady-state covariance matrix
#' for the latent variables.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param theta Numeric matrix.
#'   The covariance matrix of
#'   the measurement error
#'   (\eqn{\boldsymbol{\Theta}}).
#' @param cov_eta Numeric matrix.
#'   The steady-state covariance matrix
#'   for the latent variables
#'   in the state space model
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
SSMCovY <- function(lambda, theta, cov_eta) {
  .SSMCovY(
    lambda = lambda,
    theta = theta,
    cov_eta = cov_eta
  )
}
