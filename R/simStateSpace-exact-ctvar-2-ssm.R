#' Convert CTVAR Parameters to State Space Parameters
#'
#' This function converts continuous-time vector autoregressive (CTVAR)
#' parameters to discrete-time state space model parameters by exact
#' sampling of the derivative companion state.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param iota Numeric vector of length \eqn{q}.
#'   Constant term entering the highest derivative of the CTVAR model.
#' @param phi List of length \eqn{p}.
#'   Numeric \eqn{q \times q} matrices containing the CTVAR derivative
#'   coefficient matrices \eqn{\mathbf{S}_{0}, \ldots,
#'   \mathbf{S}_{p - 1}}, in increasing derivative order.
#' @param sigma_l Numeric \eqn{q \times q} matrix.
#'   Lower triangular factor of the continuous-time diffusion covariance:
#'   `sigma_l %*% t(sigma_l)`. Noise enters the highest derivative.
#'   For a positive-definite covariance matrix `sigma`, use `t(chol(sigma))`.
#'   Zero and positive-semidefinite covariances are allowed.
#' @param delta_t Positive finite numeric scalar.
#'   Sampling interval \eqn{\Delta t}.
#'
#' @details
#' The transformation uses Equations (4)--(7) and Table 1, block A,
#' of Chow et al. (2023). For \eqn{q} observed variables and derivative
#' order \eqn{p}, the state has length \eqn{pq} and is ordered as
#' \eqn{(\mathbf{y}, \mathbf{y}^{(1)}, \ldots, \mathbf{y}^{(p-1)})}.
#' The drift coefficients are \eqn{\mathbf{S}_0, \ldots, \mathbf{S}_{p-1}};
#' the intercept and diffusion enter the last block of the state.
#'
#' For orders above one, the sampled observed process generally requires
#' a state space representation. For order one,
#' this transformation corresponds to
#' [LinSDE2SSM()] with `phi = phi[[1]]` and an identity observation matrix.
#'
#' @return Returns a list of state space model parameters:
#'   - `alpha`: Numeric vector of length \eqn{pq}.
#'     Intercept of the sampled derivative state.
#'   - `beta`: Numeric \eqn{pq \times pq} matrix.
#'     Transition matrix of the sampled derivative state.
#'   - `psi_l`: Numeric \eqn{pq \times pq} matrix.
#'     Lower triangular factor of the interval process noise covariance,
#'     reconstructed as `psi_l %*% t(psi_l)`.
#'   - `lambda`: Numeric \eqn{q \times pq} matrix.
#'     Observation matrix selecting the first \eqn{q} state components.
#'
#' @references
#'   Chow, S.-M., Losardo, D., Park, J. J., & Molenaar, P. C. M. (2023).
#'   Continuous-time dynamic models: Connections to structural equation models
#'   and other discrete-time models.
#'   In R. H. Hoyle (Ed.),
#'   _Handbook of structural equation modeling_ (2nd ed.).
#'   The Guilford Press.
#'
#' @examples
#' # Example 1: CTVAR(1) to a state space model
#' iota <- c(0.32, 0.23)
#' phi <- list(
#'   matrix(
#'     data = c(
#'       -0.10, 0.05,
#'       0.05, -0.10
#'     ),
#'     nrow = 2
#'   )
#' )
#' sigma <- matrix(
#'   data = c(
#'     2.79, 0.06,
#'     0.06, 3.27
#'   ),
#'   nrow = 2
#' )
#' sigma_l <- t(chol(sigma))
#' CTVAR2SSM(
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   delta_t = 1.0
#' )
#'
#' # Example 2: CTVAR(2) to a state space model
#' iota <- c(0.32, 0.23)
#' phi <- list(
#'   matrix(
#'     data = c(
#'       -0.10, 0.05,
#'       0.05, -0.10
#'     ),
#'     nrow = 2
#'   ),
#'   matrix(
#'     data = c(
#'       -0.30, 0.02,
#'       0.02, -0.30
#'     ),
#'     nrow = 2
#'   )
#' )
#' sigma <- matrix(
#'   data = c(
#'     2.79, 0.06,
#'     0.06, 3.27
#'   ),
#'   nrow = 2
#' )
#' sigma_l <- t(chol(sigma))
#' CTVAR2SSM(
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   delta_t = 1.0
#' )
#'
#' @seealso [LinSDE2SSM()], [SSM2LinSDE()]
#' @family Parameter Transformation Functions
#' @keywords simStateSpace transformation ctvar ssm
#' @export
CTVAR2SSM <- function(iota,
                      phi,
                      sigma_l,
                      delta_t = 1) {
  .CheckCTVAR(
    iota = iota,
    phi = phi,
    sigma_l = sigma_l,
    delta_t = delta_t
  )
  .CTVAR2SSM(
    iota = c(iota),
    phi = phi,
    sigma_l = sigma_l,
    delta_t = delta_t
  )
}
