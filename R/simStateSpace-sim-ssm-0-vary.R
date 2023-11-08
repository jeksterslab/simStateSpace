#' Simulate Data using a State Space Model Parameterization
#' for n > 1 Individuals (Varying Parameters)
#'
#' This function simulates data
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
#'   `alpha`,
#'   `beta`,
#'   `psi_sqrt`,
#'   `nu`,
#'   `lambda`, and
#'   `theta_sqrt`)
#'   is less the `n`,
#'   the function will cycled through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu0 List of numeric vectors.
#'   Mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#' @param sigma0_sqrt List of numeric matrices.
#'   Cholesky decomposition of the covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#' @param alpha List of numeric vectors.
#'   Vector of intercepts for the dynamic model
#'   (\eqn{\boldsymbol{\alpha}}).
#' @param beta List of numeric matrices.
#'   Transition matrix relating the values of the latent variables
#'   at time `t - 1` to those at time `t`
#'   (\eqn{\boldsymbol{\beta}}).
#' @param psi_sqrt List of numeric matrices.
#'   Cholesky decomposition of the process noise covariance matrix
#'   (\eqn{\boldsymbol{\Psi}}).
#' @param nu List of numeric vectors.
#'   Vector of intercepts for the measurement model
#'   (\eqn{\boldsymbol{\nu}}).
#' @param lambda List of numeric matrices.
#'   Factor loading matrix linking the latent variables
#'   to the observed variables
#'   (\eqn{\boldsymbol{\Lambda}}).
#' @param theta_sqrt List of numeric matrices.
#'   Cholesky decomposition of the measurement error covariance matrix
#'   (\eqn{\boldsymbol{\Theta}}).
#' @inheritParams SimSSM0Fixed
#' @inherit SimSSM0Fixed return
#' @inherit SimSSM0 references
#'
#' @examples
#' # prepare parameters
#' # In this example, beta varies across individuals
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- list(null_vec)
#' sigma0_sqrt <- list(iden_sqrt)
#' alpha <- list(null_vec)
#' beta <- list(
#'   diag(x = 0.1, nrow = k),
#'   diag(x = 0.2, nrow = k),
#'   diag(x = 0.3, nrow = k),
#'   diag(x = 0.4, nrow = k),
#'   diag(x = 0.5, nrow = k)
#' )
#' psi_sqrt <- list(iden_sqrt)
#' nu <- list(null_vec)
#' lambda <- list(iden)
#' theta_sqrt <- list(chol(diag(x = 0.50, nrow = k)))
#' time <- 50
#' burn_in <- 0
#'
#' ssm <- SimSSM0Vary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim
#' @export
SimSSM0Vary <- function(n,
                        mu0,
                        sigma0_sqrt,
                        alpha,
                        beta,
                        psi_sqrt,
                        nu,
                        lambda,
                        theta_sqrt,
                        time,
                        burn_in) {
  stopifnot(
    is.list(mu0),
    is.list(sigma0_sqrt),
    is.list(alpha),
    is.list(beta),
    is.list(psi_sqrt),
    is.list(nu),
    is.list(lambda),
    is.list(theta_sqrt)
  )
  stopifnot(
    length(mu0) <= n,
    length(sigma0_sqrt) <= n,
    length(alpha) <= n,
    length(beta) <= n,
    length(psi_sqrt) <= n,
    length(nu) <= n,
    length(lambda) <= n,
    length(theta_sqrt) <= n
  )
  foo <- function(i,
                  mu0,
                  sigma0_sqrt,
                  alpha,
                  beta,
                  psi_sqrt,
                  nu,
                  lambda,
                  theta_sqrt) {
    dat_i <- SimSSM0(
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
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
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      SIMPLIFY = FALSE
    )
  )
}
