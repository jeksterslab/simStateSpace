#' Intercept from
#' Steady-State Mean Vector for the
#' Latent Variables in the
#' Linear Stochastic Differential Equation Model
#'
#' The intercept vector
#' for the latent variables
#' in the linear stochastic differential equation model
#' \eqn{\boldsymbol{\iota}}
#' is given by
#' \deqn{
#'   \boldsymbol{\iota}
#'   =
#'   - \boldsymbol{\Phi} \mathrm{Mean} \left( \boldsymbol{\eta} \right)
#' }
#' where
#' \eqn{\boldsymbol{\Phi}}
#' is the drift matrix,
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is the steady-state mean vector
#' for the latent variables.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams LinSDE2SSM
#' @inheritParams LinSDEInterceptY
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
#' LinSDEMeanEta(
#'   phi = phi,
#'   iota = iota
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace linsde
#' @export
LinSDEInterceptEta <- function(phi, mean_eta) {
  .LinSDEInterceptEta(
    phi = phi,
    mean_eta = mean_eta
  )
}
