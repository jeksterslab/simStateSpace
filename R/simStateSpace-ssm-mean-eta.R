#' Steady-State Mean Vector for the
#' Latent Variables in the
#' State Space Model
#'
#' The steady-state mean vector
#' for the latent variables
#' in the state space model
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is given by
#' \deqn{
#'   \mathrm{Mean} \left( \boldsymbol{\eta} \right)
#'   =
#'   \left(
#'     \mathbf{I} - \boldsymbol{\beta}
#'   \right)^{-1}
#'   \boldsymbol{\alpha}
#' }
#' where
#' \eqn{\boldsymbol{\beta}}
#' is the transition matrix relating the values of the latent variables
#' at the previous to the current time point,
#' \eqn{\boldsymbol{\alpha}}
#' is a vector of constant values for the dynamic model,
#' and
#' \eqn{\mathbf{I}}
#' is an identity matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMFixed
#'
#' @examples
#' beta <- 0.50 * diag(3)
#' alpha <- rep(x = 0.001, times = 3)
#' SSMMeanEta(
#'   beta = beta,
#'   alpha = alpha
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SSMMeanEta <- function(beta, alpha) {
  .SSMMeanEta(
    beta = beta,
    alpha = alpha
  )
}
