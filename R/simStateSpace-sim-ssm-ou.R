#' Simulate Data from the Ornstein–Uhlenbeck Model
#' using a State Space Model Parameterization (n = 1)
#'
#' This function simulates data from the Ornstein–Uhlenbeck model
#' using a state space model parameterization.
#' See details for more information.
#'
#' @details
#'   ## Type 0
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Theta}
#'     \right)
#'   }
#'   where \eqn{\mathbf{y}_{t}}, \eqn{\boldsymbol{\eta}_{t}},
#'   and \eqn{\boldsymbol{\varepsilon}_{t}}
#'   are random variables and \eqn{\boldsymbol{\nu}},
#'   \eqn{\boldsymbol{\Lambda}},
#'   and \eqn{\boldsymbol{\Theta}} are model parameters.
#'   \eqn{\mathbf{y}_{t}} is a vector of observed random variables
#'   at time \eqn{t},
#'   \eqn{\boldsymbol{\eta}_{t}} is a vector of latent random variables
#'   at time \eqn{t},
#'   and \eqn{\boldsymbol{\varepsilon}_{t}}
#'   is a vector of random measurement errors
#'   at time \eqn{t},
#'   while \eqn{\boldsymbol{\nu}} is a vector of intercept,
#'   \eqn{\boldsymbol{\Lambda}} is a matrix of factor loadings,
#'   and \eqn{\boldsymbol{\Theta}} is the covariance matrix of
#'   \eqn{\boldsymbol{\varepsilon}}.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{t}
#'     =
#'     \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{t}
#'   }
#'   where \eqn{\boldsymbol{\mu}} is the long-term mean or equilibrium level,
#'   \eqn{\boldsymbol{\Phi}} is the rate of mean reversion,
#'   determining how quickly the variable returns to its mean,
#'   \eqn{\boldsymbol{\Sigma}} is the matrix of volatility
#'   or randomness in the process, and \eqn{\mathrm{d}\boldsymbol{W}}
#'   is a Wiener process or Brownian motion,
#'   which represents random fluctuations.
#'
#'   ## Type 1
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Theta}
#'     \right) .
#'   }
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{t}
#'     =
#'     \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{t}
#'   }
#'   where
#'   \eqn{\mathbf{x}_{t}} is a vector of covariates at time \eqn{t},
#'   and \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}} is the coefficient matrix
#'   linking the covariates to the latent variables.
#'
#'   ## Type 2
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\Gamma}_{\mathbf{y}}
#'     \mathbf{x}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Theta}
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\Gamma}_{\mathbf{y}}} is the coefficient matrix
#'   linking the covariates to the observed variables.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{t}
#'     =
#'     \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{t} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu Numeric vector.
#'   The long-term mean or equilibrium level
#'   (\eqn{\boldsymbol{\mu}}).
#' @param phi Numeric matrix.
#'   The rate of mean reversion,
#'   determining how quickly the variable returns to its mean
#'   (\eqn{\boldsymbol{\Phi}}).
#' @param sigma_sqrt Numeric matrix.
#'   Cholesky decomposition of the matrix of volatility
#'   or randomness in the process
#'   (\eqn{\boldsymbol{\Sigma}}).
#' @param delta_t Numeric.
#'   Time interval (\eqn{\delta_t}).
#' @inheritParams SimSSM
#'
#' @references
#'   Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
#'   Continuous-time dynamic models: Connections to structural equation models
#'   and other discrete-time models.
#'   In R. H. Hoyle (Ed.), Handbook of structural equation modeling (2nd ed.).
#'   The Guilford Press.
#'
#'   Uhlenbeck, G. E., & Ornstein, L. S. (1930).
#'   On the theory of the brownian motion.
#'   *Physical Review*, *36*(5), 823–841.
#'   \doi{10.1103/physrev.36.823}
#'
#' @return Returns a list with the following elements:
#'   - `y`: A `t` by `k` matrix of values for the manifest variables.
#'   - `eta`: A `t` by `p` matrix of values for the latent variables.
#'   - `time`: A vector of continuous time points of length `t`
#'      starting from 0 with `delta_t` increments.
#'   - `n`: Number of individuals.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' p <- k <- 2
#' iden <- diag(p)
#' iden_sqrt <- chol(iden)
#' mu0 <- c(-3.0, 1.5)
#' sigma0_sqrt <- iden_sqrt
#' mu <- c(5.76, 5.18)
#' phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
#' sigma_sqrt <- chol(
#'   matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
#' )
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta_sqrt <- chol(diag(x = 0.50, nrow = k))
#' delta_t <- 0.10
#' time <- 50
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
#' x <- matrix(
#'   data = rnorm(n = k * (time + burn_in)),
#'   ncol = k
#' )
#'
#' # Type 0
#' ssm <- SimSSMOU(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   type = 0,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 1
#' ssm <- SimSSMOU(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 2
#' ssm <- SimSSMOU(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ou
#' @export
SimSSMOU <- function(mu0,
                     sigma0_sqrt,
                     mu,
                     phi,
                     sigma_sqrt,
                     nu,
                     lambda,
                     theta_sqrt,
                     gamma_y = NULL,
                     gamma_eta = NULL,
                     x = NULL,
                     type = 0,
                     delta_t,
                     time,
                     burn_in) {
  stopifnot(
    type %in% 0:2
  )
  if (type == 0) {
    return(
      .SimSSM0OU(
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        mu = mu,
        phi = phi,
        sigma_sqrt = sigma_sqrt,
        nu = nu,
        lambda = lambda,
        theta_sqrt = theta_sqrt,
        delta_t = delta_t,
        time = time,
        burn_in = burn_in
      )
    )
  }
  if (type == 1) {
    return(
      .SimSSM1OU(
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        mu = mu,
        phi = phi,
        sigma_sqrt = sigma_sqrt,
        nu = nu,
        lambda = lambda,
        theta_sqrt = theta_sqrt,
        gamma_eta = gamma_eta,
        x = x,
        delta_t = delta_t,
        time = time,
        burn_in = burn_in
      )
    )
  }
  if (type == 2) {
    return(
      .SimSSM2OU(
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        mu = mu,
        phi = phi,
        sigma_sqrt = sigma_sqrt,
        nu = nu,
        lambda = lambda,
        theta_sqrt = theta_sqrt,
        gamma_y = gamma_y,
        gamma_eta = gamma_eta,
        x = x,
        delta_t = delta_t,
        time = time,
        burn_in = burn_in
      )
    )
  }
}
