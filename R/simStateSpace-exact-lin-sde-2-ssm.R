#' Convert Linear SDE Parameters to State Space Parameters
#'
#' This function converts linear stochastic differential equation (SDE)
#' parameters to discrete-time state space model parameters by exact
#' discretization over a specified sampling interval.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param iota Numeric vector of length \eqn{q}.
#'   An unobserved term that is constant over time
#'   (\eqn{\boldsymbol{\iota}}).
#' @param phi Numeric \eqn{q \times q} matrix.
#'   The drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations
#'   (\eqn{\boldsymbol{\Phi}}).
#' @param sigma_l Numeric \eqn{q \times q} matrix.
#'   Lower triangular factor of the continuous-time diffusion covariance:
#'   `sigma_l %*% t(sigma_l)`. For a positive-definite covariance matrix
#'   `sigma`, use `t(chol(sigma))`.
#'   The factor must be lower triangular with a non-negative diagonal.
#'   Zero and positive-semidefinite covariances are allowed.
#' @param delta_t Positive finite numeric scalar.
#'   Sampling interval \eqn{\Delta t}.
#'
#' @details
#' The linear stochastic differential equation is
#' \deqn{
#'   \mathrm{d}\boldsymbol{\eta}(t)
#'   = (\boldsymbol{\iota} + \boldsymbol{\Phi}\boldsymbol{\eta}(t))
#'     \mathrm{d}t + \mathbf{L}\,\mathrm{d}\mathbf{W}(t),
#' }
#' where \eqn{\mathbf{L}} is `sigma_l` and
#' \eqn{\boldsymbol{\Sigma} = \mathbf{L}\mathbf{L}^{\prime}}.
#' Over a sampling interval \eqn{\Delta t}, the state equation is
#' \deqn{
#'   \boldsymbol{\eta}(t + \Delta t)
#'   = \boldsymbol{\alpha} + \boldsymbol{\beta}\boldsymbol{\eta}(t)
#'     + \boldsymbol{\zeta}(t),
#'   \qquad \boldsymbol{\zeta}(t) \sim
#'   \mathcal{N}(\mathbf{0}, \boldsymbol{\Psi}).
#' }
#' The exact discrete-time parameters are
#' \deqn{
#'   \boldsymbol{\beta} = \exp(\Delta t\boldsymbol{\Phi}),
#'   \qquad
#'   \boldsymbol{\alpha} = \left[\int_0^{\Delta t}
#'     \exp(s\boldsymbol{\Phi})\,\mathrm{d}s\right]\boldsymbol{\iota},
#' }
#' and
#' \deqn{
#'   \boldsymbol{\Psi} = \int_0^{\Delta t}
#'     \exp(s\boldsymbol{\Phi})\boldsymbol{\Sigma}
#'     \exp(s\boldsymbol{\Phi}^{\prime})\,\mathrm{d}s.
#' }
#' The implementation evaluates these integrals using block matrix
#' exponentials, so zero and singular drift matrices are supported.
#' All inputs must contain finite values and have compatible dimensions.
#'
#' @return Returns a list of state space parameters:
#'   - `alpha`: Numeric vector.
#'     Vector of constant values for the dynamic model
#'     (\eqn{\boldsymbol{\alpha}}).
#'   - `beta`: Numeric matrix.
#'     Transition matrix relating the values of the latent variables
#'     from the previous time point to the current time point
#'     (\eqn{\boldsymbol{\beta}}).
#'   - `psi_l`: Numeric matrix.
#'     Lower triangular factor
#'     of the process noise covariance matrix
#'     \eqn{\boldsymbol{\Psi}}.
#'
#' @references
#'   Harvey, A. C. (1990).
#'   Forecasting, structural time series models and the Kalman filter.
#'   Cambridge University Press.
#'   \doi{10.1017/cbo9781107049994}
#'
#' @examples
#' # Example 1: Exact sampling of a bivariate linear SDE
#' iota <- c(0.32, 0.23)
#'
#' phi <- matrix(
#'   data = c(
#'     -0.10, 0.05,
#'     0.05, -0.10
#'   ),
#'   nrow = 2
#' )
#'
#' sigma <- matrix(
#'   data = c(
#'     2.79, 0.06,
#'     0.06, 3.27
#'   ),
#'   nrow = 2
#' )
#'
#' sigma_l <- t(chol(sigma))
#' delta_t <- 1.0
#'
#' LinSDE2SSM(
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   delta_t = delta_t
#' )
#'
#' # Example 2: Zero drift (Brownian motion with a constant drift term)
#' LinSDE2SSM(
#'   iota = 0.2,
#'   phi = matrix(0),
#'   sigma_l = matrix(0.5),
#'   delta_t = 0.5
#' )
#'
#' @seealso [SSM2LinSDE()], [CTVAR2SSM()]
#' @family Parameter Transformation Functions
#' @keywords simStateSpace transformation linsde
#' @export
LinSDE2SSM <- function(iota,
                       phi,
                       sigma_l,
                       delta_t = 1.0) {
  tol <- sqrt(.Machine$double.eps)
  if (is.matrix(iota)) {
    dim(iota) <- NULL
  }
  stopifnot(
    "`phi` must be a matrix." = is.matrix(phi),
    "`phi` must be numeric." = is.numeric(phi),
    "`phi` must contain only finite values." = all(is.finite(phi)),
    "`phi` must have positive dimensions." = nrow(phi) > 0L,
    "`phi` must be a square matrix." = nrow(phi) == ncol(phi),
    "`iota` must be a vector." = is.null(dim(iota)),
    "`iota` must be numeric." = is.numeric(iota),
    "`iota` must contain only finite values." = all(is.finite(iota)),
    "`nrow(phi)` and `length(iota)` must match." = identical(
      x = nrow(phi),
      y = length(iota)
    ),
    "`sigma_l` must be a matrix." = is.matrix(sigma_l),
    "`sigma_l` must be numeric." = is.numeric(sigma_l),
    "`sigma_l` must contain only finite values." = all(is.finite(sigma_l)),
    "`delta_t` must be numeric." = is.numeric(delta_t),
    "`delta_t` must have length 1." = length(delta_t) == 1L,
    "`delta_t` must be finite." = is.finite(delta_t),
    "`delta_t` must be positive." = delta_t > 0.0
  )
  q <- nrow(phi)
  stopifnot(
    "`sigma_l` must have dimensions compatible with `phi`." =
      identical(
        x = dim(sigma_l),
        y = c(q, q)
      ),
    "`sigma_l` must be lower triangular." = all(
      abs(sigma_l[upper.tri(sigma_l)]) <= tol
    ),
    "`sigma_l` must have non-negative diagonal values." = all(
      diag(sigma_l) >= -tol
    )
  )
  .LinSDE2SSM(
    iota = iota,
    phi = phi,
    sigma_l = sigma_l,
    delta_t = delta_t
  )
}
