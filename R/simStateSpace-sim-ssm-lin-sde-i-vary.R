#' Simulate Data from the
#' Linear Stochastic Differential Equation Model
#' using a State Space Model Parameterization
#' (Individual-Varying Parameters)
#'
#' This function simulates data from the
#' linear stochastic differential equation model
#' using a state space model parameterization.
#' It assumes that the parameters can vary
#' across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0_l`,
#'   `iota`,
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
#' @param type Integer.
#'   State space model type.
#'   See Details in [SimSSMLinSDEFixed()] for more information.
#' @inheritParams SimSSMIVary
#' @param iota List of numeric vectors.
#'   Each element of the list
#'   is an unobserved term that is constant over time
#'   (\eqn{\boldsymbol{\iota}}).
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
#' @inherit SimSSMFixed return
#' @inherit SimSSMLinSDEFixed references
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
#' sigma0 <- 0.001 * diag(p)
#' sigma0_l <- list(
#'   t(chol(sigma0))
#' )
#' iota <- list(
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
#' theta <- 0.001 * diag(k)
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
#' ssm <- SimSSMLinSDEIVary(
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
#' ssm <- SimSSMLinSDEIVary(
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
#' ssm <- SimSSMLinSDEIVary(
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
SimSSMLinSDEIVary <- function(n, time, delta_t = 1.0,
                              mu0, sigma0_l,
                              iota, phi, sigma_l,
                              nu, lambda, theta_l,
                              type = 0,
                              x = NULL, gamma = NULL, kappa = NULL) {
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  call <- match.call()
  args <- list(
    n = n, time = time, delta_t = delta_t,
    mu0 = mu0, sigma0_l = sigma0_l,
    iota = iota, phi = phi, sigma_l = sigma_l,
    nu = nu, lambda = lambda, theta_l = theta_l,
    type = type,
    x = x, gamma = gamma, kappa = kappa,
    ou = FALSE
  )
  model <- list(
    model = "linsde",
    covariates = covariates,
    fixed = FALSE,
    vary_i = TRUE
  )
  fun <- "SimSSMLinSDEIVary"
  .CheckSimulation(
    n = n, time = time, delta_t = delta_t,
    par = list(
      mu0 = mu0, sigma0_l = sigma0_l,
      iota = iota, phi = phi, sigma_l = sigma_l,
      nu = nu, lambda = lambda, theta_l = theta_l
    ),
    type = type, types = 0:2,
    x = x, gamma = gamma, kappa = kappa,
    fixed = FALSE
  )
  ssm <- .LinSDE2SSMSets(
    n = n,
    iota = iota,
    phi = phi,
    sigma_l = sigma_l,
    delta_t = delta_t,
    ou = FALSE
  )
  data <- .SimSSMData(
    n = n,
    time = time,
    delta_t = delta_t,
    mu0 = mu0, sigma0_l = sigma0_l,
    alpha = ssm$alpha, beta = ssm$beta, psi_l = ssm$psi_l,
    nu = nu, lambda = lambda, theta_l = theta_l,
    type = type,
    x = x, gamma = gamma, kappa = kappa,
    fixed = FALSE
  )
  out <- list(
    call = call,
    args = args,
    model = model,
    data = data,
    fun = fun
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  out
}
