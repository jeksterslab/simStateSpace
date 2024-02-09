#' Simulate Data from the
#' Ornstein–Uhlenbeck Model
#' using a State Space Model Parameterization
#' (Individual-Varying Parameters)
#'
#' This function simulates data from the
#' Ornstein–Uhlenbeck model
#' using a state space model parameterization.
#' In this model,
#' the parameters can vary across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0_l`,
#'   `mu`,
#'   `phi`,
#'   `sigma_l`,
#'   `nu`,
#'   `lambda`,
#'   `theta_l`,
#'   `gamma`, or
#'   `kappa`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMIVary
#' @param mu List of numeric vectors.
#'   Each element of the list
#'   is the long-term mean or equilibrium level
#'   (\eqn{\boldsymbol{\mu}}).
#' @param phi List of numeric matrix.
#'   Each element of the list
#'   is the drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations
#'   (\eqn{\boldsymbol{\Phi}}).
#'   The negative value of `phi` is the rate of mean reversion,
#'   determining how quickly the variable returns to its mean
#'   (\eqn{- \boldsymbol{\Phi}}).
#' @param sigma_l List of numeric matrix.
#'   Each element of the list
#'   is the Cholesky factorization (`t(chol(sigma))`)
#'   of the covariance matrix of volatility
#'   or randomness in the process
#'   \eqn{\boldsymbol{\Sigma}}.
#'
#' @inherit SimSSMFixed references return
#' @inheritParams SimSSMLinSDEIVary
#'
#' @examples
#' # prepare parameters
#' # In this example, phi varies across individuals.
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' delta_t <- 0.10
#' ## dynamic structure
#' p <- 2
#' mu0 <- list(
#'   c(-3.0, 1.5)
#' )
#' sigma0 <- diag(p)
#' sigma0_l <- list(
#'   t(chol(sigma0))
#' )
#' mu <- list(
#'   c(5.76, 5.18)
#' )
#' phi <- list(
#'   -0.1 * diag(p),
#'   -0.2 * diag(p),
#'   -0.3 * diag(p),
#'   -0.4 * diag(p),
#'   -0.5 * diag(p)
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
#' sigma_l <- list(
#'   t(chol(sigma))
#' )
#' ## measurement model
#' k <- 2
#' nu <- list(
#'   rep(x = 0, times = k)
#' )
#' lambda <- list(
#'   diag(k)
#' )
#' theta <- 0.50 * diag(k)
#' theta_l <- list(
#'   t(chol(theta))
#' )
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
#' gamma <- list(
#'   diag(x = 0.10, nrow = p, ncol = j)
#' )
#' kappa <- list(
#'   diag(x = 0.10, nrow = k, ncol = j)
#' )
#'
#' # Type 0
#' ssm <- SimSSMOUIVary(
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
#' ssm <- SimSSMOUIVary(
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
#' ssm <- SimSSMOUIVary(
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
SimSSMOUIVary <- function(n, time, delta_t = 1.0,
                          mu0, sigma0_l,
                          mu, phi, sigma_l,
                          nu, lambda, theta_l,
                          type = 0,
                          x = NULL, gamma = NULL, kappa = NULL) {
  stopifnot(type %in% c(0, 1, 2))
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  if (type == 0) {
    data <- .SimSSMLinSDEIVary0(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      iota = rep(x = mu, length.out = n),
      phi = rep(x = phi, length.out = n),
      sigma_l = rep(x = sigma_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      ou = TRUE
    )
  }
  if (type == 1) {
    stopifnot(
      !is.null(x),
      !is.null(gamma)
    )
    data <- .SimSSMLinSDEIVary1(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      iota = rep(x = mu, length.out = n),
      phi = rep(x = phi, length.out = n),
      sigma_l = rep(x = sigma_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      x = rep(x = x, length.out = n),
      gamma = rep(x = gamma, length.out = n),
      ou = TRUE
    )
  }
  if (type == 2) {
    stopifnot(
      !is.null(x),
      !is.null(gamma),
      !is.null(kappa)
    )
    data <- .SimSSMLinSDEIVary2(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      iota = rep(x = mu, length.out = n),
      phi = rep(x = phi, length.out = n),
      sigma_l = rep(x = sigma_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      x = rep(x = x, length.out = n),
      gamma = rep(x = gamma, length.out = n),
      kappa = rep(x = kappa, length.out = n),
      ou = TRUE
    )
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n, time = time,
      mu0 = mu0, sigma0_l = sigma0_l,
      iota = mu, phi = phi, sigma_l = sigma_l,
      nu = nu, lambda = lambda, theta_l = theta_l,
      type = type,
      x = x, gamma = gamma, kappa = kappa,
      ou = TRUE
    ),
    model = list(
      model = "ou",
      covariates = covariates,
      fixed = FALSE,
      vary_i = TRUE
    ),
    data = data,
    fun = "SimSSMOUIVary"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
