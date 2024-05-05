#' Simulate Data from the
#' Ornstein<U+2013>Uhlenbeck Model
#' using a State Space Model Parameterization
#' (Fixed Parameters)
#'
#' This function simulates data from the
#' Ornstein<U+2013>Uhlenbeck (OU) model
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
#'     - \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t}
#'   }
#'   where
#'   \eqn{\boldsymbol{\mu}}
#'   is the long-term mean or equilibrium level,
#'   \eqn{- \boldsymbol{\Phi}}
#'   is the rate of mean reversion,
#'   determining how quickly the variable returns to its mean,
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
#'     - \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
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
#'     - \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
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
#' ## The OU model as a linear stochastic differential equation model
#'
#'   The OU model is a first-order
#'   linear stochastic differential equation model
#'   in the form of
#'
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
#'   \eqn{\boldsymbol{\mu} = - \boldsymbol{\Phi}^{-1} \boldsymbol{\iota}}
#'   and, equivalently
#'   \eqn{\boldsymbol{\iota} = - \boldsymbol{\Phi} \boldsymbol{\mu}}.
#'
#' @references
#'   Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
#'   Equivalence and differences between structural equation modeling
#'   and state-space modeling techniques.
#'   *Structural Equation Modeling: A Multidisciplinary Journal*,
#'   17(2), 303<U+2013>332.
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
#'   Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011).
#'   A hierarchical latent stochastic differential equation model
#'   for affective dynamics.
#'   Psychological Methods,
#'   16 (4), 468<U+2013>490.
#'   \doi{10.1037/a0024375}
#'
#'   Uhlenbeck, G. E., & Ornstein, L. S. (1930).
#'   On the theory of the brownian motion.
#'   Physical Review,
#'   36 (5), 823<U+2013>841.
#'   \doi{10.1103/physrev.36.823}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu Numeric vector.
#'   The long-term mean or equilibrium level
#'   (\eqn{\boldsymbol{\mu}}).
#' @param phi Numeric matrix.
#'   The drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations
#'   (\eqn{\boldsymbol{\Phi}}).
#'   The negative value of `phi` is the rate of mean reversion,
#'   determining how quickly the variable returns to its mean
#'   (\eqn{- \boldsymbol{\Phi}}).
#' @inheritParams SimSSMLinSDEFixed
#' @inherit SimSSMFixed return
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
#' mu <- c(5.76, 5.18)
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
#' ssm <- SimSSMOUFixed(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   mu = mu,
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
#' ssm <- SimSSMOUFixed(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   mu = mu,
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
#' ssm <- SimSSMOUFixed(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   mu = mu,
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
#' @keywords simStateSpace sim ou
#' @export
SimSSMOUFixed <- function(n, time, delta_t = 1.0,
                          mu0, sigma0_l,
                          mu, phi, sigma_l,
                          nu, lambda, theta_l,
                          type = 0,
                          x = NULL, gamma = NULL, kappa = NULL) {
  iota <- -phi %*% mu
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
    mu = mu,
    iota = iota,
    phi = phi,
    sigma_l = sigma_l
  )
  out$model$model <- "ou"
  out$fun <- "SimSSMOUFixed"
  return(out)
}
