#' Simulate Data from the State Space Model
#' (Fixed Parameters)
#'
#' This function simulates data from the
#' state space model.
#' In this model,
#' the parameters are invariant cross individuals and across time.
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
#'   is a vector of observed random variables,
#'   \eqn{\boldsymbol{\eta}_{i, t}}
#'   is a vector of latent random variables,
#'   and
#'   \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   is a vector of random measurement errors,
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\boldsymbol{\nu}}
#'   is a vector of intercepts,
#'   \eqn{\boldsymbol{\Lambda}}
#'   is a matrix of factor loadings,
#'   and
#'   \eqn{\boldsymbol{\Theta}}
#'   is the covariance matrix of
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
#'     \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{i, t - 1}
#'     +
#'     \boldsymbol{\zeta}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\eta}_{i, t}},
#'   \eqn{\boldsymbol{\eta}_{i, t - 1}},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{i, t}}
#'   are random variables,
#'   and
#'   \eqn{\boldsymbol{\alpha}},
#'   \eqn{\boldsymbol{\beta}},
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   are model parameters.
#'   \eqn{\boldsymbol{\eta}_{i, t}}
#'   is a vector of latent variables
#'   at time \eqn{t} and individual \eqn{i},
#'   \eqn{\boldsymbol{\eta}_{i, t - 1}}
#'   is a vector of latent variables
#'   at time \eqn{t - 1} and individual \eqn{i},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{i, t}}
#'   is a vector of dynamic noise
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\boldsymbol{\alpha}}
#'   is a vector of intercepts,
#'   \eqn{\boldsymbol{\beta}}
#'   is a matrix of autoregression
#'   and cross regression coefficients,
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   is the covariance matrix of
#'   \eqn{\boldsymbol{\zeta}_{i, t}}.
#'
#'   An alternative representation of the dynamic noise
#'   is given by
#'   \deqn{
#'     \boldsymbol{\zeta}_{i, t}
#'     =
#'     \boldsymbol{\Psi}^{\frac{1}{2}}
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
#'   \eqn{
#'     \left( \boldsymbol{\Psi}^{\frac{1}{2}} \right)
#'     \left( \boldsymbol{\Psi}^{\frac{1}{2}} \right)^{\prime}
#'     =
#'     \boldsymbol{\Psi} .
#'   }
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
#'     \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{i, t - 1}
#'     +
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\zeta}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right)
#'   }
#'   where
#'   \eqn{\mathbf{x}_{i, t}} is a vector of covariates
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\Gamma}} is the coefficient matrix
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
#'     \boldsymbol{\Kappa}
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
#'   \eqn{\boldsymbol{\Kappa}} is the coefficient matrix
#'   linking the covariates to the observed variables.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{i, t - 1}
#'     +
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\zeta}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right) .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of individuals.
#' @param time Positive integer.
#'   Number of time points.
#' @param delta_t Numeric.
#'   Time interval.
#'   The default value is `1.0`
#'   with an option to use a numeric value
#'   for the discretized state space model
#'   parameterization of the
#'   linear stochastic differential equation model.
#' @param mu0 Numeric vector.
#'   Mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#' @param sigma0_l Numeric matrix.
#'   Cholesky factorization (`t(chol(sigma0))`)
#'   of the covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#' @param alpha Numeric vector.
#'   Vector of constant values for the dynamic model
#'   (\eqn{\boldsymbol{\alpha}}).
#' @param beta Numeric matrix.
#'   Transition matrix relating the values of the latent variables
#'   at the previous to the current time point
#'   (\eqn{\boldsymbol{\beta}}).
#' @param psi_l Numeric matrix.
#'   Cholesky factorization (`t(chol(psi))`)
#'   of the covariance matrix
#'   of the process noise
#'   (\eqn{\boldsymbol{\Psi}}).
#' @param nu Numeric vector.
#'   Vector of intercept values for the measurement model
#'   (\eqn{\boldsymbol{\nu}}).
#' @param lambda Numeric matrix.
#'   Factor loading matrix linking the latent variables
#'   to the observed variables
#'   (\eqn{\boldsymbol{\Lambda}}).
#' @param theta_l Numeric matrix.
#'   Cholesky factorization (`t(chol(theta))`)
#'   of the covariance matrix
#'   of the measurement error
#'   (\eqn{\boldsymbol{\Theta}}).
#' @param type Integer.
#'   State space model type.
#'   See Details for more information.
#' @param x List.
#'   Each element of the list is a matrix of covariates
#'   for each individual `i` in `n`.
#'   The number of columns in each matrix
#'   should be equal to `time`.
#' @param gamma Numeric matrix.
#'   Matrix linking the covariates to the latent variables
#'   at current time point
#'   (\eqn{\boldsymbol{\Gamma}}).
#' @param kappa Numeric matrix.
#'   Matrix linking the covariates to the observed variables
#'   at current time point
#'   (\eqn{\boldsymbol{\Kappa}}).
#'
#' @references
#'   Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
#'   Equivalence and differences between structural equation modeling
#'   and state-space modeling techniques.
#'   *Structural Equation Modeling: A Multidisciplinary Journal*,
#'   17(2), 303–332.
#'   \doi{10.1080/10705511003661553}
#'
#' @return Returns an object of class `simstatespace`
#'   which is a list with the following elements:
#'   - `call`: Function call.
#'   - `args`: Function arguments.
#'   - `data`: Generated data which is a list of length `n`.
#'     Each element of `data` is a list with the following elements:
#'     * `id`: A vector of ID numbers with length `t`,
#'       where `t` is the value of the function argument `time`.
#'     * `time`: A vector time points of length `t`.
#'     * `y`: A `t` by `k` matrix of values for the manifest variables.
#'     * `eta`: A `t` by `p` matrix of values for the latent variables.
#'     * `x`: A `t` by `j` matrix of values for the covariates
#'       (when covariates are included).
#'   - `fun`: Function used.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- diag(p)
#' psi_l <- t(chol(psi))
#' ## measurement model
#' k <- 3
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.50 * diag(k)
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
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
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
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
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
#' @keywords simStateSpace sim ssm
#' @export
SimSSMFixed <- function(n, time, delta_t = 1.0,
                        mu0, sigma0_l,
                        alpha, beta, psi_l,
                        nu, lambda, theta_l,
                        type = 0,
                        x = NULL, gamma = NULL, kappa = NULL) {
  stopifnot(type %in% c(0, 1, 2))
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  if (type == 0) {
    data <- .SimSSMFixed0(
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      nu = nu, lambda = lambda, theta_l = theta_l
    )
  }
  if (type == 1) {
    stopifnot(
      !is.null(x),
      !is.null(gamma)
    )
    data <- .SimSSMFixed1(
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      nu = nu, lambda = lambda, theta_l = theta_l,
      x = x, gamma = gamma
    )
  }
  if (type == 2) {
    stopifnot(
      !is.null(x),
      !is.null(gamma),
      !is.null(kappa)
    )
    data <- .SimSSMFixed2(
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      nu = nu, lambda = lambda, theta_l = theta_l,
      x = x, gamma = gamma, kappa = kappa
    )
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n, time = time,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      nu = nu, lambda = lambda, theta_l = theta_l,
      type = type,
      x = x, gamma = gamma, kappa = kappa
    ),
    model = list(
      model = "ssm",
      covariates = covariates,
      fixed = TRUE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMFixed"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
