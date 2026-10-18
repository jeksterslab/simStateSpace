#' Steady-State Mean Vector for the
#' Latent Variables in the
#' Linear Stochastic Differential Equation Model
#'
#' The steady-state mean vector
#' for the latent variables
#' in the linear stochastic differential equation model
#' \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
#' is given by
#' \deqn{
#'   \mathrm{Mean} \left( \boldsymbol{\eta} \right)
#'   =
#'   -\boldsymbol{\Phi}^{-1} \boldsymbol{\iota}
#' }
#' where
#' \eqn{\boldsymbol{\Phi}}
#' is the drift matrix,
#' and
#' \eqn{\boldsymbol{\iota}}
#' is an unobserved term that is constant over time.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams LinSDE2SSM
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
#' LinSDEMeanEta(
#'   phi = phi,
#'   iota = iota
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace linsde
#' @export
LinSDEMeanEta <- function(phi,
                          iota) {
  if (is.matrix(iota)) {
    dim(iota) <- NULL
  }
  stopifnot(
    "`phi` must be a matrix." = is.matrix(phi),
    "`phi` must be numeric." = is.numeric(phi),
    "`phi` must contain only finite values." = all(is.finite(phi)),
    "`phi` must be a square matrix." = nrow(phi) == ncol(phi),
    "`iota` must be a vector." = is.null(dim(iota)),
    "`iota` must be numeric." = is.numeric(iota),
    "`iota` must contain only finite values." = all(is.finite(iota)),
    "`nrow(phi)` and `length(iota)` must match." = identical(
      x = nrow(phi),
      y = length(iota)
    )
  )
  .LinSDEMeanEta(
    phi = phi,
    iota = iota
  )
}
