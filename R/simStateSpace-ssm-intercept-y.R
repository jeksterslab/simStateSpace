#' Intercept from
#' Steady-State Mean Vector for the
#' Observed Variables in the
#' State Space Model
#'
#' The intercept vector
#' for the observed variables
#' in the state space model
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
#' @param mean_y Numeric vector.
#'   Steady-state mean vector
#'   of the observed variables
#'   \eqn{\mathrm{Mean} \left( \mathbf{y} \right)}.
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
#' mean_y <- SSMMeanY(
#'   nu = nu,
#'   lambda = lambda,
#'   mean_eta = mean_eta
#' )
#' SSMInterceptY(
#'   mean_y = mean_y,
#'   mean_eta = mean_eta,
#'   lambda = lambda
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SSMInterceptY <- function(mean_y, mean_eta, lambda) {
  .SSMInterceptY(
    mean_y = mean_y,
    mean_eta = mean_eta,
    lambda = lambda
  )
}
