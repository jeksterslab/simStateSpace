#' Steady-State Covariance Matrix for the
#' Linear Stochastic Differential Equation Model
#'
#' The steady-state covariance matrix is the solution
#' to the Sylvester equation,
#' i.e.
#' \deqn{
#'   \mathbf{A} \mathbf{X} +
#'   \mathbf{X} \mathbf{B} +
#'   \mathbf{C} = \mathbf{0} ,
#' } where \eqn{\mathbf{X}} is unknown,
#' \eqn{\mathbf{A} = \boldsymbol{\Phi}},
#' \eqn{\mathbf{B} = \boldsymbol{\Phi}^{\prime}}, and
#' \eqn{\mathbf{C} = \boldsymbol{\Sigma}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param sigma Numeric matrix.
#'   The covariance matrix of volatility
#'   or randomness in the process
#'   (\eqn{\boldsymbol{\Sigma}}).
#' @inheritParams LinSDE2SSM
#'
#' @examples
#' phi <- matrix(
#'   data = c(
#'     -0.10,
#'     0.05,
#'     0.05,
#'     -0.10
#'   ),
#'   nrow = 2
#' )
#' sigma <- matrix(
#'   data = c(
#'     2.79,
#'     0.06,
#'     0.06,
#'     3.27
#'   ),
#'   nrow = 2
#' )
#' LinSDECov(phi = phi, sigma = sigma)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim linsde
#' @export
LinSDECov <- function(phi, sigma) {
  .SolveSyl(
    A = phi,
    B = t(phi),
    C = sigma
  )
}
