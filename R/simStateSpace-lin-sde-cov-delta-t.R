#' Interval Covariance Matrix
#' in the Linear Stochastic Differential Equation Model
#'
#' The covariance matrix of the process noise
#' accumulated over a sampling interval
#' \eqn{\Delta t}
#' in the linear stochastic differential equation model,
#' i.e.
#' \deqn{
#'   \boldsymbol{\Psi}_{\Delta t} =
#'   \int_{0}^{\Delta t}
#'   \exp \left( s \boldsymbol{\Phi} \right)
#'   \boldsymbol{\Sigma}
#'   \exp \left( s \boldsymbol{\Phi}^{\prime} \right)
#'   \mathrm{d} s ,
#' }
#' where
#' \eqn{\boldsymbol{\Phi}}
#' is the drift matrix
#' and
#' \eqn{\boldsymbol{\Sigma}}
#' is the covariance matrix
#' of volatility or randomness in the process.
#'
#' @details
#' The integral is computed with the block matrix exponential
#' of Van Loan (1978),
#' composed over subintervals
#' so that strongly decaying drifts do not overflow.
#' The drift may be singular
#' (for example, \eqn{\boldsymbol{\Phi} = \mathbf{0}}
#' gives \eqn{\Delta t \boldsymbol{\Sigma}}).
#' The map is linear in \eqn{\boldsymbol{\Sigma}},
#' so `sigma` may be any symmetric matrix,
#' not only a covariance matrix.
#' For a positive semidefinite `sigma`,
#' the result is the process noise covariance
#' `psi_l %*% t(psi_l)` returned by [LinSDE2SSM()].
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param sigma Numeric matrix.
#'   Symmetric matrix,
#'   usually the covariance matrix of volatility
#'   or randomness in the process
#'   (\eqn{\boldsymbol{\Sigma}}).
#' @inheritParams LinSDE2SSM
#'
#' @return Returns a symmetric numeric matrix
#'   (\eqn{\boldsymbol{\Psi}_{\Delta t}}).
#'
#' @references
#'   Van Loan, C. F. (1978).
#'   Computing integrals involving the matrix exponential.
#'   *IEEE Transactions on Automatic Control*, *23*(3), 395--404.
#'   \doi{10.1109/TAC.1978.1101743}
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
#' sigma <- matrix(
#'   data = c(
#'     0.24455556, 0.02201587, -0.05004762,
#'     0.02201587, 0.07067800, 0.01539456,
#'     -0.05004762, 0.01539456, 0.07553061
#'   ),
#'   nrow = 3
#' )
#' LinSDECovDeltaT(
#'   phi = phi,
#'   sigma = sigma,
#'   delta_t = 1.0
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace linsde
#' @export
LinSDECovDeltaT <- function(phi,
                            sigma,
                            delta_t = 1.0) {
  stopifnot(
    "`phi` must be a matrix." = is.matrix(phi),
    "`phi` must be numeric." = is.numeric(phi),
    "`phi` must contain only finite values." = all(is.finite(phi)),
    "`phi` must have positive dimensions." = nrow(phi) > 0L,
    "`phi` must be a square matrix." = nrow(phi) == ncol(phi),
    "`sigma` must be a matrix." = is.matrix(sigma),
    "`sigma` must be numeric." = is.numeric(sigma),
    "`sigma` must contain only finite values." = all(is.finite(sigma)),
    "`phi` and `sigma` must have the same dimensions." = identical(
      x = dim(phi),
      y = dim(sigma)
    ),
    "`sigma` must be symmetric." = isTRUE(
      all.equal(
        target = sigma,
        current = t(sigma),
        tolerance = sqrt(.Machine$double.eps)
      )
    ),
    "`delta_t` must be numeric." = is.numeric(delta_t),
    "`delta_t` must have length 1." = length(delta_t) == 1L,
    "`delta_t` must be finite." = is.finite(delta_t),
    "`delta_t` must be positive." = delta_t > 0.0
  )
  .LinSDECovDeltaT(
    phi = phi,
    sigma = (sigma + t(sigma)) / 2,
    delta_t = delta_t
  )
}
