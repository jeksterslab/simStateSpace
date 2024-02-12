#' Simulate Data from the
#' Linear Growth Curve Model
#'
#' This function simulates data from the
#' linear growth curve model.
#'
#' @details
#'   ## Type 0
#'
#'   The measurement model is given by
#'   \deqn{
#'     Y_{i, t}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 0 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'       \eta_{0_{i, t}} \\
#'       \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta
#'     \right)
#'   }
#'   where \eqn{Y_{i, t}}, \eqn{\eta_{0_{i, t}}},
#'   \eqn{\eta_{1_{i, t}}},
#'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   are random variables and
#'   \eqn{\theta} is a model parameter.
#'   \eqn{Y_{i, t}} is the observed random variable
#'   at time \eqn{t} and individual \eqn{i},
#'   \eqn{\eta_{0_{i, t}}} (intercept)
#'   and
#'   \eqn{\eta_{1_{i, t}}} (slope)
#'   form a vector of latent random variables
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   a vector of random measurement errors
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\theta} is the variance of
#'   \eqn{\boldsymbol{\varepsilon}}.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 1 \\
#'       0 & 1 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t - 1}} \\
#'     \eta_{1_{i, t - 1}} \\
#'     \end{array}
#'     \right) .
#'   }
#'
#'   The mean vector and covariance matrix of the intercept and slope
#'   are captured in the mean vector and covariance matrix
#'   of the initial condition given by
#'   \deqn{
#'     \boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}
#'     =
#'     \left(
#'     \begin{array}{c}
#'       \mu_{\eta_{0}} \\
#'       \mu_{\eta_{1}} \\
#'     \end{array}
#'     \right) \quad \mathrm{and,}
#'   }
#'
#'   \deqn{
#'     \boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       \sigma^{2}_{\eta_{0}} &
#'       \sigma_{\eta_{0}, \eta_{1}} \\
#'       \sigma_{\eta_{1}, \eta_{0}} &
#'       \sigma^{2}_{\eta_{1}} \\
#'     \end{array}
#'     \right) .
#'   }
#'
#'   ## Type 1
#'
#'   The measurement model is given by
#'   \deqn{
#'     Y_{i, t}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 0 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'       \eta_{0_{i, t}} \\
#'       \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta
#'     \right) .
#'   }
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 1 \\
#'       0 & 1 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t - 1}} \\
#'     \eta_{1_{i, t - 1}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t}
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
#'     Y_{i, t}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 0 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'       \eta_{0_{i, t}} \\
#'       \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
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
#'     0,
#'     \theta
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\kappa}} represents the coefficient matrix
#'   linking the covariates to the observed variables.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 1 \\
#'       0 & 1 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t - 1}} \\
#'     \eta_{1_{i, t - 1}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu0 Numeric vector.
#'   A vector of length two.
#'   The first element is the mean of the intercept,
#'   and the second element is the mean of the slope.
#' @param sigma0_l Numeric matrix.
#'   Cholesky factorization (`t(chol(sigma0))`)
#'   of the covariance matrix
#'   of the intercept and the slope.
#' @param theta_l Numeric.
#'   Square root of the common measurement error variance.
#'
#' @inheritParams SimSSMFixed
#'
#' @inherit SimSSMFixed references return
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 5
#' ## dynamic structure
#' p <- 2
#' mu0 <- c(0.615, 1.006)
#' sigma0 <- matrix(
#'   data = c(
#'     1.932,
#'     0.618,
#'     0.618,
#'     0.587
#'   ),
#'   nrow = p
#' )
#' sigma0_l <- t(chol(sigma0))
#' ## measurement model
#' k <- 1
#' theta <- 0.50
#' theta_l <- sqrt(theta)
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = j * time),
#'         nrow = j
#'       )
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#' kappa <- diag(x = 0.10, nrow = k, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
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
#' @keywords simStateSpace sim growth
#' @export
SimSSMLinGrowth <- function(n, time,
                            mu0, sigma0_l, theta_l,
                            type = 0,
                            x = NULL, gamma = NULL, kappa = NULL) {
  stopifnot(type %in% c(0, 1, 2))
  p <- 2
  k <- 1
  stopifnot(
    length(mu0) == p,
    dim(sigma0_l) == c(p, p)
  )
  theta_l <- as.matrix(
    theta_l
  )
  stopifnot(
    dim(theta_l) == c(k, k)
  )
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  alpha <- rep(x = 0, times = p)
  beta <- matrix(
    data = c(1, 0, 1, 1),
    nrow = p
  )
  psi_l <- matrix(
    data = 0,
    nrow = p,
    ncol = p
  )
  nu <- rep(x = 0, times = k)
  lambda <- matrix(
    data = c(1, 0),
    nrow = k
  )
  if (type == 0) {
    data <- .SimSSMFixed0(
      n = n,
      time = time,
      delta_t = 1.0,
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
      delta_t = 1.0,
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
      delta_t = 1.0,
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
      model = "lingrowth",
      covariates = covariates,
      fixed = TRUE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMLinGrowth"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
