#' Steady-State Mean Vector for the
#' Observed Variables in the
#' Linear Stochastic Differential Equation Model
#'
#' The steady-state mean vector
#' for the observed variables
#' in the linear stochastic differential equation model
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is given by
#' \deqn{
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
#' iota <- c(0.317, 0.230)
#' phi <- matrix(
#'   data = c(
#'     -0.10, 0.05,
#'     0.05, -0.10
#'   ),
#'   nrow = 2
#' )
#' nu <- rep(x = 0.03, times = 2)
#' lambda <- diag(2)
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
