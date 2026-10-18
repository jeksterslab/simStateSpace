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
#' is the drift matrix and
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
LinSDEInterceptEta <- function(phi,
                               mean_eta) {
  if (is.matrix(mean_eta)) {
    dim(mean_eta) <- NULL
  }
  stopifnot(
    "`phi` must be a matrix." = is.matrix(phi),
    "`phi` must be numeric." = is.numeric(phi),
    "`phi` must contain only finite values." = all(is.finite(phi)),
    "`phi` must be a square matrix." = nrow(phi) == ncol(phi),
    "`mean_eta` must be a vector." = is.null(dim(mean_eta)),
    "`mean_eta` must be numeric." = is.numeric(mean_eta),
    "`mean_eta` must contain only finite values." = all(is.finite(mean_eta)),
    "`nrow(phi)` and `length(mean_eta)` must match." = identical(
      x = nrow(phi),
      y = length(mean_eta)
    )
  )
  .LinSDEInterceptEta(
    phi = phi,
    mean_eta = mean_eta
  )
}
