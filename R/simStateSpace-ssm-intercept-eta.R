#' Intercept from
#' Steady-State Mean Vector for the
#' Latent Variables in the
#' State Space Model
#'
#' The intercept vector
#' for the latent variables
#' in the state space model
#' \eqn{\boldsymbol{\alpha}}
#' is given by
#' \deqn{
#'   \boldsymbol{\alpha}
#'   =
#'   \mathrm{Mean} \left( \boldsymbol{\eta} \right) -
#'   \boldsymbol{\beta} \mathrm{Mean} \left( \boldsymbol{\eta} \right)
#' }
#' where
#' \eqn{\boldsymbol{\beta}}
#' is the transition matrix relating the values of the latent variables
#' at the previous to the current time point,
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is the steady-state mean vector
#' for the latent variables.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
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
#' mean_eta <- SSMMeanEta(
#'   beta = beta,
#'   alpha = alpha
#' )
#' SSMInterceptEta(
#'   beta = beta,
#'   mean_eta = mean_eta
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SSMInterceptEta <- function(beta, mean_eta) {
  .SSMInterceptEta(
    beta = beta,
    mean_eta = mean_eta
  )
}
