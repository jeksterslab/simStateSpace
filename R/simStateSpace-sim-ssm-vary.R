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
#'   `lambda`,
#'   `theta_sqrt`,
#'   `gamma_y`, or
#'   `gamma_eta`)
#'   is less the `n`,
#'   the function will cycle through the available values.
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
#' @inheritParams SimSSMFixed
#' @inherit SimSSMFixed return
#' @inherit SimSSM references
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
#' gamma_y <- gamma_eta <- list(0.10 * diag(k))
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMVary(
#'   n = n,
#'   type = 0,
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
#' # Type 1
#' ssm <- SimSSMVary(
#'   n = n,
#'   type = 1,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 2
#' ssm <- SimSSMVary(
#'   n = n,
#'   type = 2,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ssm
#' @export
SimSSMVary <- function(n,
                       type,
                       mu0,
                       sigma0_sqrt,
                       alpha,
                       beta,
                       psi_sqrt,
                       nu,
                       lambda,
                       theta_sqrt,
                       gamma_y = NULL,
                       gamma_eta = NULL,
                       x = NULL,
                       time = 0,
                       burn_in = 0) {
  stopifnot(
    type %in% 0:2
  )
  if (type == 0) {
    return(
      .SimSSM0Vary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_sqrt = rep(x = psi_sqrt, length.out = n),
        nu = rep(x = nu, length.out = n),
        lambda = rep(x = lambda, length.out = n),
        theta_sqrt = rep(x = theta_sqrt, length.out = n),
        time = time,
        burn_in = burn_in
      )
    )
  }
  if (type == 1) {
    return(
      .SimSSM1Vary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_sqrt = rep(x = psi_sqrt, length.out = n),
        nu = rep(x = nu, length.out = n),
        lambda = rep(x = lambda, length.out = n),
        theta_sqrt = rep(x = theta_sqrt, length.out = n),
        gamma_eta = rep(x = gamma_eta, length.out = n),
        x = x,
        time = time,
        burn_in = burn_in
      )
    )
  }
  if (type == 2) {
    return(
      .SimSSM2Vary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_sqrt = rep(x = psi_sqrt, length.out = n),
        nu = rep(x = nu, length.out = n),
        lambda = rep(x = lambda, length.out = n),
        theta_sqrt = rep(x = theta_sqrt, length.out = n),
        gamma_y = rep(x = gamma_y, length.out = n),
        gamma_eta = rep(x = gamma_eta, length.out = n),
        x = x,
        time = time,
        burn_in = burn_in
      )
    )
  }
}
