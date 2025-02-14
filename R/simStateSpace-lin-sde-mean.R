#' Steady-State Mean Vector for the
#' Linear Stochastic Differential Equation Model
#'
#' The steady-state mean vector is given by
#' \deqn{
#'   -\boldsymbol{\Phi}^{-1} \boldsymbol{\iota}
#' }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams LinSDE2SSM
#'
#' @examples
#' iota <- c(0.317, 0.230)
#' phi <- matrix(
#'   data = c(
#'     -0.10,
#'     0.05,
#'     0.05,
#'     -0.10
#'   ),
#'   nrow = 2
#' )
#' LinSDEMean(phi = phi, iota = iota)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim linsde
#' @export
LinSDEMean <- function(phi, iota) {
  solve(-phi, iota)
}
