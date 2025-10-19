#' Steady-State Mean Vector for the
#' Observed Variables in the
#' State Space Model
#'
#' The steady-state mean vector
#' for the observed variables
#' in the state space model
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
#' beta <- matrix(
#'   data = c(
#'     0.7, 0.5, -0.1,
#'     0.0, 0.6, 0.4,
#'     0.0, 0.0, 0.5
#'   ),
#'   nrow = 3
#' )
#' alpha <- rep(x = 1, times = 3)
#' lambda <- diag(3)
#' nu <- rep(x = 1, times = 3)
#' mean_eta <- SSMMeanEta(
#'   beta = beta,
#'   alpha = alpha
#' )
#' SSMMeanY(
#'   nu = nu,
#'   lambda = lambda,
#'   mean_eta = mean_eta
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SSMMeanY <- function(nu, lambda, mean_eta) {
  .SSMMeanY(
    nu = nu,
    lambda = lambda,
    mean_eta = mean_eta
  )
}
