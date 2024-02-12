#' Simulate Data from the
#' Linear Stochastic Differential Equation Model
#' using a State Space Model Parameterization
#' (Fixed Parameters)
#'
#' This function simulates data from the
#' linear stochastic differential equation model
#' using a state space model parameterization.
#' It assumes that the parameters remain constant
#' across individuals and over time.
#'
#' @details
#'   ## Type 0
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{i, t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{i, t}
#'     +
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Theta}
#'     \right)
#'   }
#'   where
#'   \eqn{\mathbf{y}_{i, t}},
#'   \eqn{\boldsymbol{\eta}_{i, t}},
#'   and
#'   \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   are random variables
#'   and
#'   \eqn{\boldsymbol{\nu}},
#'   \eqn{\boldsymbol{\Lambda}},
#'   and
#'   \eqn{\boldsymbol{\Theta}}
#'   are model parameters.
#'   \eqn{\mathbf{y}_{i, t}}
#'   represents a vector of observed random variables,
#'   \eqn{\boldsymbol{\eta}_{i, t}}
#'   a vector of latent random variables,
#'   and
#'   \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   a vector of random measurement errors,
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\boldsymbol{\nu}}
#'   denotes a vector of intercepts,
#'   \eqn{\boldsymbol{\Lambda}}
#'   a matrix of factor loadings,
#'   and
#'   \eqn{\boldsymbol{\Theta}}
#'   the covariance matrix of
#'   \eqn{\boldsymbol{\varepsilon}}.
#'
#'   An alternative representation of the measurement error
#'   is given by
#'   \deqn{
#'     \boldsymbol{\varepsilon}_{i, t}
#'     =
#'     \boldsymbol{\Theta}^{\frac{1}{2}}
#'     \mathbf{z}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \mathbf{z}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \mathbf{I}
#'     \right)
#'   }
#'   where
#'   \eqn{\mathbf{z}_{i, t}} is a vector of
#'   independent standard normal random variables and
#'   \eqn{
#'     \left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)
#'     \left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)^{\prime}
#'     =
#'     \boldsymbol{\Theta} .
#'   }
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{i, t}
#'     =
#'     \left(
#'     \boldsymbol{\iota}
#'     +
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t}
#'   }
#'   where
#'   \eqn{\boldsymbol{\iota}}
#'   is a term which is unobserved and constant over time,
#'   \eqn{\boldsymbol{\Phi}}
#'   is the drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations,
#'   \eqn{\boldsymbol{\Sigma}}
#'   is the matrix of volatility
#'   or randomness in the process, and
#'   \eqn{\mathrm{d}\boldsymbol{W}}
#'   is a Wiener process or Brownian motion,
#'   which represents random fluctuations.
#'
#'   ## Type 1
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{i, t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{i, t}
#'     +
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
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
#'     \mathrm{d} \boldsymbol{\eta}_{i, t}
#'     =
#'     \left(
#'     \boldsymbol{\iota}
#'     +
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t}
#'   }
#'   where
#'   \eqn{\mathbf{x}_{i, t}} represents a vector of covariates
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\Gamma}} the coefficient matrix
#'   linking the covariates to the latent variables.
#'
#'   ## Type 2
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{i, t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{i, t}
#'     +
#'     \boldsymbol{\kappa}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Theta}
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\kappa}} represents the coefficient matrix
#'   linking the covariates to the observed variables.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{i, t}
#'     =
#'     \left(
#'     \boldsymbol{\iota}
#'     +
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t} .
#'   }
#'
#' ## State Space Parameterization
#'
#'   The state space parameters
#'   as a function of the linear stochastic differential equation model
#'   parameters
#'   are given by
#'   \deqn{
#'       \boldsymbol{\beta}
#'       =
#'       \exp{
#'         \left(
#'           \boldsymbol{\Phi}
#'           \Delta_{t}
#'         \right)
#'       }
#'   }
#'
#'   \deqn{
#'       \boldsymbol{\alpha}
#'       =
#'       \boldsymbol{\Phi}^{-1}
#'       \left(
#'         \boldsymbol{\beta} - \mathbf{I}_{p}
#'       \right)
#'       \boldsymbol{\iota}
#'   }
#'
#'   \deqn{
#'       \mathrm{vec}
#'       \left(
#'         \boldsymbol{\Psi}
#'       \right)
#'       =
#'       \left[
#'         \left(
#'           \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
#'         \right)
#'         +
#'         \left(
#'           \mathbf{I}_{p} \otimes \boldsymbol{\Phi}
#'         \right)
#'       \right]
#'       \left[
#'         \exp
#'         \left(
#'           \left[
#'             \left(
#'               \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
#'             \right)
#'             +
#'             \left(
#'               \mathbf{I}_{p} \otimes \boldsymbol{\Phi}
#'             \right)
#'           \right]
#'           \Delta_{t}
#'         \right)
#'         -
#'         \mathbf{I}_{p \times p}
#'       \right]
#'       \mathrm{vec}
#'       \left(
#'         \boldsymbol{\Sigma}
#'       \right)
#'   }
#'   where \eqn{p} is the number of latent variables and
#'   \eqn{\Delta_{t}} is the time interval.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams LinSDE2SSM
#' @inheritParams SimSSMFixed
#' @inherit SimSSMFixed references return
#'
#' @references
#'   Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
#'   Equivalence and differences between structural equation modeling
#'   and state-space modeling techniques.
#'   *Structural Equation Modeling: A Multidisciplinary Journal*,
#'   17(2), 303–332.
#'   \doi{10.1080/10705511003661553}
#'
#'   Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
#'   Continuous-time dynamic models:
#'   Connections to structural equation models and other discrete-time models.
#'   In R. H. Hoyle (Ed.),
#'   Handbook of structural equation modeling (2nd ed.).
#'   The Guilford Press.
#'
#'   Harvey, A. C. (1990).
#'   Forecasting, structural time series models and the Kalman filter.
#'   Cambridge University Press.
#'   \doi{10.1017/cbo9781107049994}
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' delta_t <- 0.10
#' ## dynamic structure
#' p <- 2
#' mu0 <- c(-3.0, 1.5)
#' sigma0 <- 0.001 * diag(p)
#' sigma0_l <- t(chol(sigma0))
#' iota <- c(0.317, 0.230)
#' phi <- matrix(
#'   data = c(
#'     -0.10,
#'     0.05,
#'     0.05,
#'     -0.10
#'   ),
#'   nrow = p
#' )
#' sigma <- matrix(
#'   data = c(
#'     2.79,
#'     0.06,
#'     0.06,
#'     3.27
#'   ),
#'   nrow = p
#' )
#' sigma_l <- t(chol(sigma))
#' ## measurement model
#' k <- 2
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.001 * diag(k)
#' theta_l <- t(chol(theta))
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     matrix(
#'       data = stats::rnorm(n = time * j),
#'       nrow = j,
#'       ncol = time
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#' kappa <- diag(x = 0.10, nrow = k, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMLinSDEFixed(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMLinSDEFixed(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinSDEFixed(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma = gamma,
#'   kappa = kappa
#' )
#'
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim linsde
#' @export
SimSSMLinSDEFixed <- function(n, time, delta_t = 1.0,
                              mu0, sigma0_l,
                              iota, phi, sigma_l,
                              nu, lambda, theta_l,
                              type = 0,
                              x = NULL, gamma = NULL, kappa = NULL) {
  ssm <- LinSDE2SSM(
    iota = iota,
    phi = phi,
    sigma_l = sigma_l,
    delta_t = delta_t
  )
  out <- SimSSMFixed(
    n = n, time = time, delta_t = delta_t,
    mu0 = mu0, sigma0_l = sigma0_l,
    alpha = ssm$alpha, beta = ssm$beta, psi_l = ssm$psi_l,
    nu = nu, lambda = lambda, theta_l = theta_l,
    type = type,
    x = x, gamma = gamma, kappa = kappa
  )
  out$args <- c(
    out$args,
    iota = iota,
    phi = phi,
    sigma_l = sigma_l
  )
  out$model$model <- "linsde"
  out$fun <- "SimSSMLinSDEFixed"
  return(out)
}
