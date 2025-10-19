#' Steady-State Mean Vector for the
#' Observed Variables in the
#' Linear Stochastic Differential Equation Model
#'
#' The steady-state mean vector
#' for the observed variables
#' in the linear stochastic differential equation model
#' \eqn{\mathrm{Mean} \left( \mathbf{y} \right)}
#' is given by
#' \deqn{
#'   \mathrm{Mean} \left( \mathbf{y} \right)
#'   =
#'   \boldsymbol{\nu}
#'   +
#'   \boldsymbol{\Lambda}
#'   \mathrm{Mean} \left( \boldsymbol{\eta} \right)
#' }
#' where
#' \eqn{\boldsymbol{\nu}}
#' is the vector of intercept values
#' for the measurement model,
#' \eqn{\boldsymbol{\Lambda}}
#' is the matrix of factor loadings,
#' and
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is the steady-state mean vector
#' for the latent variables.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mean_eta Numeric vector.
#'   Steady-state mean vector
#'   of the latent variables
#'   \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}.
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
#' iota <- rep(x = 1, times = 3)
#' lambda <- diag(3)
#' nu <- rep(x = 1, times = 3)
#' mean_eta <- LinSDEMeanEta(
#'   phi = phi,
#'   iota = iota
#' )
#' LinSDEMeanY(
#'   nu = nu,
#'   lambda = lambda,
#'   mean_eta = mean_eta
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace linsde
#' @export
LinSDEMeanY <- function(nu, lambda, mean_eta) {
  .LinSDEMeanY(
    nu = nu,
    lambda = lambda,
    mean_eta = mean_eta
  )
}
