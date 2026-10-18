#' Convert State Space Parameters to Linear SDE Parameters
#'
#' This function converts parameters from
#' the discrete-time state space model parameterization
#' to the linear stochastic differential equation model
#' using the exact matrix-logarithm mapping.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param alpha Numeric vector of length \eqn{q}.
#'   Intercept of the discrete-time state equation.
#' @param beta Numeric \eqn{q \times q} matrix.
#'   Transition matrix of the discrete-time state equation.
#' @param psi_l Numeric \eqn{q \times q} matrix.
#'   Lower triangular factor of the interval process noise covariance:
#'   `psi_l %*% t(psi_l)`. For a positive-definite covariance matrix `psi`,
#'   use `t(chol(psi))`. Zero and positive-semidefinite covariances are allowed.
#' @param delta_t Positive finite numeric scalar.
#'   Sampling interval \eqn{\Delta t}.
#' @param tol Non-negative finite numeric scalar.
#'   Tolerance for the largest absolute imaginary component of the
#'   matrix logarithm.
#' @param tol_psi_l Non-negative finite numeric scalar.
#'   Tolerance for checking that `psi_l` is lower triangular
#'   with a non-negative diagonal.
#'
#' @details
#' The state equation and its continuous-time counterpart are described in
#' [LinSDE2SSM()]. The drift matrix is computed using the principal matrix
#' logarithm:
#' \deqn{
#'   \boldsymbol{\Phi} = \frac{1}{\Delta t}\log(\boldsymbol{\beta}).
#' }
#' The intercept and diffusion covariance are recovered by solving
#' \deqn{
#'   \boldsymbol{\alpha} = \left[\int_0^{\Delta t}
#'     \exp(s\boldsymbol{\Phi})\,\mathrm{d}s\right]\boldsymbol{\iota}
#' }
#' and
#' \deqn{
#'   \boldsymbol{\Psi} = \int_0^{\Delta t}
#'     \exp(s\boldsymbol{\Phi})\boldsymbol{\Sigma}
#'     \exp(s\boldsymbol{\Phi}^{\prime})\,\mathrm{d}s,
#' }
#' where `psi_l %*% t(psi_l)` gives \eqn{\boldsymbol{\Psi}}.
#' The integrated exponential operators support zero and singular drifts.
#'
#' The principal matrix-logarithm branch does not search all possible
#' continuous-time embeddings or resolve temporal aliasing. Inputs are
#' rejected if the logarithm is not finite (singular `beta`),
#' has an imaginary component exceeding `tol`,
#' the integral inversions are not unique, or the recovered diffusion
#' covariance is not positive semidefinite.
#' All inputs must contain finite values and have compatible dimensions.
#'
#' @return Returns a list of linear stochastic differential equation
#'   model parameters:
#'   - `iota`: Numeric vector.
#'     An unobserved term that is constant over time
#'     \eqn{\boldsymbol{\iota}}.
#'   - `phi`: Numeric matrix.
#'     The drift matrix which represents the rate of change of the
#'     solution in the absence of any random fluctuations
#'     \eqn{\boldsymbol{\Phi}}.
#'   - `sigma_l`: Numeric matrix.
#'     Lower triangular factor
#'     of the covariance matrix of volatility or randomness in the process
#'     \eqn{\boldsymbol{\Sigma}}.
#'
#' @references
#'   Harvey, A. C. (1990).
#'   Forecasting, structural time series models and the Kalman filter.
#'   Cambridge University Press.
#'   \doi{10.1017/cbo9781107049994}
#'
#' @examples
#' # Example 1: Recover the parameters of a bivariate linear SDE
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
#' ssm <- LinSDE2SSM(
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   delta_t = delta_t
#' )
#'
#' SSM2LinSDE(
#'   alpha = ssm$alpha,
#'   beta = ssm$beta,
#'   psi_l = ssm$psi_l,
#'   delta_t = delta_t
#' )
#'
#' @seealso [LinSDE2SSM()], [CTVAR2SSM()]
#' @family Parameter Transformation Functions
#' @keywords simStateSpace transformation linsde ssm
#' @export
SSM2LinSDE <- function(alpha,
                       beta,
                       psi_l,
                       delta_t,
                       tol = sqrt(.Machine$double.eps),
                       tol_psi_l = sqrt(.Machine$double.eps)) {
  stopifnot(
    "`tol` must be a finite non-negative numeric scalar." =
      is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol >= 0,
    "`tol_psi_l` must be a finite non-negative numeric scalar." = all(
      is.numeric(tol_psi_l),
      length(tol_psi_l) == 1L,
      is.finite(tol_psi_l),
      tol_psi_l >= 0
    )
  )
  if (is.matrix(alpha)) {
    dim(alpha) <- NULL
  }
  stopifnot(
    "`beta` must be a matrix." = is.matrix(beta),
    "`beta` must be numeric." = is.numeric(beta),
    "`beta` must contain only finite values." = all(is.finite(beta)),
    "`beta` must have positive dimensions." = nrow(beta) > 0L,
    "`beta` must be a square matrix." = nrow(beta) == ncol(beta),
    "`alpha` must be a vector." = is.null(dim(alpha)),
    "`alpha` must be numeric." = is.numeric(alpha),
    "`alpha` must contain only finite values." = all(is.finite(alpha)),
    "`nrow(beta)` and `length(alpha)` must match." = identical(
      x = nrow(beta),
      y = length(alpha)
    ),
    "`psi_l` must be a matrix." = is.matrix(psi_l),
    "`psi_l` must be numeric." = is.numeric(psi_l),
    "`psi_l` must contain only finite values." = all(is.finite(psi_l)),
    "`delta_t` must be numeric." = is.numeric(delta_t),
    "`delta_t` must have length 1." = length(delta_t) == 1L,
    "`delta_t` must be finite." = is.finite(delta_t),
    "`delta_t` must be positive." = delta_t > 0.0
  )
  q <- nrow(beta)
  stopifnot(
    "`psi_l` must have dimensions compatible with `beta`." =
      identical(
        x = dim(psi_l),
        y = c(q, q)
      ),
    "`psi_l` must be lower triangular." = all(
      abs(psi_l[upper.tri(psi_l)]) <= tol_psi_l
    ),
    "`psi_l` must have non-negative diagonal values." = all(
      diag(psi_l) >= -tol_psi_l
    )
  )
  .SSM2LinSDE(
    alpha = alpha,
    beta = beta,
    psi_l = psi_l,
    delta_t = delta_t,
    tol = tol
  )
}
