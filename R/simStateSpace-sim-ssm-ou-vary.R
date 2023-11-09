#' Simulate Data from an Ornstein–Uhlenbeck Model
#' using a State Space Model Parameterization
#' for n > 1 Individuals (Varying Parameters)
#'
#' This function simulates data from an Ornstein–Uhlenbeck model
#' using a state space model parameterization
#' for `n > 1` individuals.
#' In this model,
#' the parameters can vary across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0_sqrt`,
#'   `mu`,
#'   `phi`,
#'   `sigma_sqrt`,
#'   `nu`,
#'   `lambda`,
#'   `theta_sqrt`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu List of numeric vectors.
#'   The long-term mean or equilibrium level
#'   (\eqn{\boldsymbol{\mu}}).
#' @param phi List of numeric matrices.
#'   The rate of mean reversion,
#'   determining how quickly the variable returns to its mean
#'   (\eqn{\boldsymbol{\Phi}}).
#' @param sigma_sqrt List of numeric matrices.
#'   Cholesky decomposition of the matrix of volatility
#'   or randomness in the process
#'   (\eqn{\boldsymbol{\Sigma}}).
#' @inheritParams SimSSMOUFixed
#' @inherit SimSSM0Fixed return
#' @inherit SimSSMOU references
#'
#' @examples
#' # prepare parameters
#' # In this example, phi varies across individuals
#' set.seed(42)
#' p <- k <- 2
#' iden <- diag(p)
#' iden_sqrt <- chol(iden)
#' n <- 5
#' mu0 <- list(c(-3.0, 1.5))
#' sigma0_sqrt <- list(iden_sqrt)
#' mu <- list(c(5.76, 5.18))
#' phi <- list(
#'   as.matrix(Matrix::expm(diag(x = -0.1, nrow = k))),
#'   as.matrix(Matrix::expm(diag(x = -0.2, nrow = k))),
#'   as.matrix(Matrix::expm(diag(x = -0.3, nrow = k))),
#'   as.matrix(Matrix::expm(diag(x = -0.4, nrow = k))),
#'   as.matrix(Matrix::expm(diag(x = -0.5, nrow = k)))
#' )
#' sigma_sqrt <- list(
#'   chol(
#'     matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
#'   )
#' )
#' nu <- list(rep(x = 0, times = k))
#' lambda <- list(diag(k))
#' theta_sqrt <- list(chol(diag(x = 0.50, nrow = k)))
#' delta_t <- 0.10
#' time <- 50
#' burn_in <- 0
#'
#' ssm <- SimSSMOUVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim
#' @export
SimSSMOUVary <- function(n,
                         mu0,
                         sigma0_sqrt,
                         mu,
                         phi,
                         sigma_sqrt,
                         nu,
                         lambda,
                         theta_sqrt,
                         delta_t,
                         time,
                         burn_in) {
  stopifnot(
    is.list(mu0),
    is.list(sigma0_sqrt),
    is.list(mu),
    is.list(phi),
    is.list(sigma_sqrt),
    is.list(nu),
    is.list(lambda),
    is.list(theta_sqrt)
  )
  stopifnot(
    length(mu0) <= n,
    length(sigma0_sqrt) <= n,
    length(mu) <= n,
    length(phi) <= n,
    length(sigma_sqrt) <= n,
    length(nu) <= n,
    length(lambda) <= n,
    length(theta_sqrt) <= n
  )
  foo <- function(i,
                  mu0,
                  sigma0_sqrt,
                  mu,
                  phi,
                  sigma_sqrt,
                  nu,
                  lambda,
                  theta_sqrt) {
    dat_i <- SimSSMOU(
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
    return(
      list(
        y = dat_i$y,
        eta = dat_i$eta,
        time = dat_i$time,
        id = matrix(data = i, ncol = 1, nrow = time)
      )
    )
  }
  return(
    mapply(
      FUN = foo,
      i = seq_len(n),
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      mu = mu,
      phi = phi,
      sigma_sqrt = sigma_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      SIMPLIFY = FALSE
    )
  )
}
