#' Simulate Data from an Ornstein–Uhlenbeck Model
#' using a State Space Model Parameterization
#' for n > 1 Individuals (Individual-Varying Parameters)
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
#'   `theta_sqrt`,
#'   `gamma_y`, or
#'   `gamma_eta`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu List of numeric vectors.
#'   Each element of the list
#'   is the long-term mean or equilibrium level
#'   (\eqn{\boldsymbol{\mu}}).
#' @param phi List of numeric matrices.
#'   Each element of the list
#'   is the rate of mean reversion,
#'   determining how quickly the variable returns to its mean
#'   (\eqn{\boldsymbol{\Phi}}).
#' @param sigma_sqrt List of numeric matrices.
#'   Each element of the list
#'   is the Cholesky decomposition of the matrix of volatility
#'   or randomness in the process
#'   (\eqn{\boldsymbol{\Sigma}}).
#'
#' @inheritParams SimSSMOUFixed
#' @inherit SimSSMFixed return
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
#' gamma_y <- gamma_eta <- list(0.10 * diag(k))
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'    )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMOUIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   type = 0,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 1
#' ssm <- SimSSMOUIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 2
#' ssm <- SimSSMOUIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ou
#' @export
SimSSMOUIVary <- function(n,
                         mu0,
                         sigma0_sqrt,
                         mu,
                         phi,
                         sigma_sqrt,
                         nu,
                         lambda,
                         theta_sqrt,
                         gamma_y = NULL,
                         gamma_eta = NULL,
                         x = NULL,
                         type = 0,
                         delta_t,
                         time,
                         burn_in = 0) {
  switch(
    EXPR = as.character(type),
    "0" = {
      return(
        .SimSSM0OUIVary(
          n = n,
          mu0 = rep(x = mu0, length.out = n),
          sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
          mu = rep(x = mu, length.out = n),
          phi = rep(x = phi, length.out = n),
          sigma_sqrt = rep(x = sigma_sqrt, length.out = n),
          nu = rep(x = nu, length.out = n),
          lambda = rep(x = lambda, length.out = n),
          theta_sqrt = rep(x = theta_sqrt, length.out = n),
          delta_t = delta_t,
          time = time,
          burn_in = burn_in
        )
      )
    },
    "1" = {
      return(
        .SimSSM1OUIVary(
          n = n,
          mu0 = rep(x = mu0, length.out = n),
          sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
          mu = rep(x = mu, length.out = n),
          phi = rep(x = phi, length.out = n),
          sigma_sqrt = rep(x = sigma_sqrt, length.out = n),
          nu = rep(x = nu, length.out = n),
          lambda = rep(x = lambda, length.out = n),
          theta_sqrt = rep(x = theta_sqrt, length.out = n),
          gamma_eta = rep(x = gamma_eta, length.out = n),
          x = x,
          delta_t = delta_t,
          time = time,
          burn_in = burn_in
        )
      )
    },
    "2" = {
      return(
        .SimSSM2OUIVary(
          n = n,
          mu0 = rep(x = mu0, length.out = n),
          sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
          mu = rep(x = mu, length.out = n),
          phi = rep(x = phi, length.out = n),
          sigma_sqrt = rep(x = sigma_sqrt, length.out = n),
          nu = rep(x = nu, length.out = n),
          lambda = rep(x = lambda, length.out = n),
          theta_sqrt = rep(x = theta_sqrt, length.out = n),
          gamma_y = rep(x = gamma_y, length.out = n),
          gamma_eta = rep(x = gamma_eta, length.out = n),
          x = x,
          delta_t = delta_t,
          time = time,
          burn_in = burn_in
        )
      )
    },
    stop(
      "Invalid `type`."
    )
  )
}
