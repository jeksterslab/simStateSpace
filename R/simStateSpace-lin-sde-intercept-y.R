#' Intercept from
#' Steady-State Mean Vector for the
#' Observed Variables in the
#' Linear Stochastic Differential Equation Model
#'
#' The intercept vector
#' for the observed variables
#' in the linear stochastic differential equation model
#' \eqn{\boldsymbol{\nu}}
#' is given by
#' \deqn{
#'   \boldsymbol{\nu}
#'   =
#'   \mathrm{Mean} \left( \mathbf{y} \right)
#'   -
#'   \boldsymbol{\Lambda}
#'   \mathrm{Mean} \left( \boldsymbol{\eta} \right)
#' }
#' where
#' \eqn{\boldsymbol{\Lambda}}
#' is the matrix of factor loadings,
#' \eqn{\mathrm{Mean} \left( \mathbf{y} \right)}
#' is the steady-state mean vector
#' for the observed variables,
#' and
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is the steady-state mean vector
#' for the latent variables.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SSMInterceptY
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
#' mean_y <- LinSDEMeanY(
#'   nu = nu,
#'   lambda = lambda,
#'   mean_eta = mean_eta
#' )
#' LinSDEInterceptY(
#'   mean_y = mean_y,
#'   mean_eta = mean_eta,
#'   lambda = lambda
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace linsde
#' @export
LinSDEInterceptY <- function(mean_y, mean_eta, lambda) {
  .SSMInterceptY(
    mean_y = mean_y,
    mean_eta = mean_eta,
    lambda = lambda
  )
}
