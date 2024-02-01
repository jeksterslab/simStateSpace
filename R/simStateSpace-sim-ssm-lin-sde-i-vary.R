#' Simulate Data from the
#' Linear Stochastic Differential Equation Model
#' using a State Space Model Parameterization
#' (Individual-Varying Parameters)
#'
#' This function simulates data from the
#' linear stochastic differential equation model
#' using a state space model parameterization.
#' In this model,
#' the parameters can vary across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0_l`,
#'   `gamma`,
#'   `phi`,
#'   `sigma_l`,
#'   `nu`,
#'   `lambda`,
#'   `theta_l`,
#'   `gamma_eta`, or
#'   `gamma_y`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMIVary
#' @param gamma List of numeric vectors.
#'   Each element of the list
#'   is an unobserved term that is constant over time
#'   (\eqn{\boldsymbol{\gamma}}).
#' @param phi List of numeric matrix.
#'   Each element of the list
#'   is the drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations
#'   (\eqn{\boldsymbol{\Phi}}).
#' @param sigma_l List of numeric matrix.
#'   Each element of the list
#'   is the Cholesky factorization (`t(chol(sigma))`)
#'   of the covariance matrix of volatility
#'   or randomness in the process
#'   \eqn{\boldsymbol{\Sigma}}.
#'
#' @inherit SimSSMFixed references return
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
#' gamma <- list(
#'   c(0.317, 0.230)
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
#' gamma_eta <- list(
#'   diag(x = 0.10, nrow = p, ncol = j)
#' )
#' gamma_y <- list(
#'   diag(x = 0.10, nrow = k, ncol = j)
#' )
#'
#' # Type 0
#' ssm <- SimSSMLinSDEIVary(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   gamma = gamma,
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
#' ssm <- SimSSMLinSDEIVary(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   gamma = gamma,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma_eta = gamma_eta
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinSDEIVary(
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   gamma = gamma,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma_eta = gamma_eta,
#'   gamma_y = gamma_y
#' )
#'
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim linsde
#' @export
SimSSMLinSDEIVary <- function(n, time, delta_t = 1.0,
                              mu0, sigma0_l,
                              gamma, phi, sigma_l,
                              nu, lambda, theta_l,
                              type = 0,
                              x = NULL, gamma_eta = NULL, gamma_y = NULL) {
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
      gamma = rep(x = gamma, length.out = n),
      phi = rep(x = phi, length.out = n),
      sigma_l = rep(x = sigma_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      ou = FALSE
    )
  }
  if (type == 1) {
    stopifnot(
      !is.null(x),
      !is.null(gamma_eta)
    )
    data <- .SimSSMLinSDEIVary1(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      gamma = rep(x = gamma, length.out = n),
      phi = rep(x = phi, length.out = n),
      sigma_l = rep(x = sigma_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      x = rep(x = x, length.out = n),
      gamma_eta = rep(x = gamma_eta, length.out = n),
      ou = FALSE
    )
  }
  if (type == 2) {
    stopifnot(
      !is.null(x),
      !is.null(gamma_eta),
      !is.null(gamma_y)
    )
    data <- .SimSSMLinSDEIVary2(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      gamma = rep(x = gamma, length.out = n),
      phi = rep(x = phi, length.out = n),
      sigma_l = rep(x = sigma_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      x = rep(x = x, length.out = n),
      gamma_eta = rep(x = gamma_eta, length.out = n),
      gamma_y = rep(x = gamma_y, length.out = n),
      ou = FALSE
    )
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n, time = time,
      mu0 = mu0, sigma0_l = sigma0_l,
      gamma = gamma, phi = phi, sigma_l = sigma_l,
      nu = nu, lambda = lambda, theta_l = theta_l,
      type = type,
      x = x, gamma_eta = gamma_eta, gamma_y = gamma_y,
      ou = FALSE
    ),
    model = list(
      model = "linsde",
      covariates = covariates,
      fixed = FALSE,
      vary_i = TRUE
    ),
    data = data,
    fun = "SimSSMLinSDEIVary"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
