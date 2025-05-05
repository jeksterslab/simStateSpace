#' Steady-State Covariance Matrix for the
#' Observed Variables in the
#' Linear Stochastic Differential Equation Model
#'
#' The steady-state covariance matrix
#' for the observed variables
#' in the linear stochastic differential equation model
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
#'   in the linear stochastic differential equation model
#' @inheritParams SimSSMFixed
#'
#' @examples
#' phi <- matrix(
#'   data = c(
#'     -0.357, 0.771, -0.450,
#'     0.0, -0.511, 0.729,
#'     0.0, 0.0, -0.693
#'   ),
#'   nrow = 3
#' )
#' sigma <- matrix(
#'   data = c(
#'     0.24455556, 0.02201587, -0.05004762,
#'     0.02201587, 0.07067800, 0.01539456,
#'     -0.05004762, 0.01539456, 0.07553061
#'   ),
#'   nrow = 3
#' )
#' lambda <- diag(3)
#' theta <- diag(3)
#' cov_eta <- LinSDECovEta(
#'   phi = phi,
#'   sigma = sigma
#' )
#' LinSDECovY(
#'   lambda = lambda,
#'   theta = theta,
#'   cov_eta = cov_eta
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace linsde
#' @export
LinSDECovY <- function(lambda, theta, cov_eta) {
  .LinSDECovY(
    lambda = lambda,
    theta = theta,
    cov_eta = cov_eta
  )
}
