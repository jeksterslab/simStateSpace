#' Simulate Data using a State Space Model Parameterization
#' for n > 1 Individuals (Individual-Varying Parameters)
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
#'   `sigma0`,
#'   `alpha`,
#'   `beta`,
#'   `psi`,
#'   `nu`,
#'   `lambda`,
#'   `theta`,
#'   `gamma_y`, or
#'   `gamma_eta`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu0 List of numeric vectors.
#'   Each element of the list
#'   is the mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#' @param sigma0 List of numeric matrices.
#'   Each element of the list
#'   is the covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#' @param alpha List of numeric vectors.
#'   Each element of the list
#'   is the vector of intercepts for the dynamic model
#'   (\eqn{\boldsymbol{\alpha}}).
#' @param beta List of numeric matrices.
#'   Each element of the list
#'   is the transition matrix relating the values of the latent variables
#'   at time `t - 1` to those at time `t`
#'   (\eqn{\boldsymbol{\beta}}).
#' @param psi List of numeric matrices.
#'   Each element of the list
#'   is the process noise covariance matrix
#'   (\eqn{\boldsymbol{\Psi}}).
#' @param nu List of numeric vectors.
#'   Each element of the list
#'   is the vector of intercepts for the measurement model
#'   (\eqn{\boldsymbol{\nu}}).
#' @param lambda List of numeric matrices.
#'   Each element of the list
#'   is the factor loading matrix linking the latent variables
#'   to the observed variables
#'   (\eqn{\boldsymbol{\Lambda}}).
#' @param theta List of numeric matrices.
#'   Each element of the list
#'   is the measurement error covariance matrix
#'   (\eqn{\boldsymbol{\Theta}}).
#'
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
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- list(null_vec)
#' sigma0 <- list(iden)
#' alpha <- list(null_vec)
#' beta <- list(
#'   diag(x = 0.1, nrow = k),
#'   diag(x = 0.2, nrow = k),
#'   diag(x = 0.3, nrow = k),
#'   diag(x = 0.4, nrow = k),
#'   diag(x = 0.5, nrow = k)
#' )
#' psi <- list(iden)
#' nu <- list(null_vec)
#' lambda <- list(iden)
#' theta <- list(diag(x = 0.50, nrow = k))
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
#' SimSSMIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' # Type 1
#' SimSSMIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' # Type 2
#' SimSSMIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ssm
#' @export
SimSSMIVary <- function(n,
                        mu0,
                        sigma0,
                        alpha,
                        beta,
                        psi,
                        nu,
                        lambda,
                        theta,
                        gamma_y = NULL,
                        gamma_eta = NULL,
                        x = NULL,
                        type,
                        time = 0,
                        burn_in = 0) {
  foo <- function(x) {
    return(
      t(chol(x))
    )
  }
  sigma0_l <- lapply(
    X = sigma0,
    FUN = foo
  )
  psi_l <- lapply(
    X = psi,
    FUN = foo
  )
  theta_l <- lapply(
    X = theta,
    FUN = foo
  )
  data <- switch(
    EXPR = as.character(type),
    "0" = {
      .SimSSM0IVary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_l = rep(x = sigma0_l, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_l = rep(x = psi_l, length.out = n),
        nu = rep(x = nu, length.out = n),
        lambda = rep(x = lambda, length.out = n),
        theta_l = rep(x = theta_l, length.out = n),
        time = time,
        burn_in = burn_in
      )
    },
    "1" = {
      .SimSSM1IVary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_l = rep(x = sigma0_l, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_l = rep(x = psi_l, length.out = n),
        nu = rep(x = nu, length.out = n),
        lambda = rep(x = lambda, length.out = n),
        theta_l = rep(x = theta_l, length.out = n),
        gamma_eta = rep(x = gamma_eta, length.out = n),
        x = x,
        time = time,
        burn_in = burn_in
      )
    },
    "2" = {
      .SimSSM2IVary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_l = rep(x = sigma0_l, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_l = rep(x = psi_l, length.out = n),
        nu = rep(x = nu, length.out = n),
        lambda = rep(x = lambda, length.out = n),
        theta_l = rep(x = theta_l, length.out = n),
        gamma_y = rep(x = gamma_y, length.out = n),
        gamma_eta = rep(x = gamma_eta, length.out = n),
        x = x,
        time = time,
        burn_in = burn_in
      )
    },
    stop(
      "Invalid `type`."
    )
  )
  if (type > 0) {
    covariates <- TRUE
  } else {
    covariates <- FALSE
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      nu = nu,
      lambda = lambda,
      theta = theta,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = type,
      time = time,
      burn_in = burn_in,
      sigma0_l = sigma0_l,
      psi_l = psi_l,
      theta_l = theta_l
    ),
    model = list(
      model = "ssm",
      n1 = FALSE,
      covariates = covariates,
      fixed = FALSE,
      vary_i = TRUE
    ),
    data = data,
    fun = "SimSSMIVary"
  )
  class(out) <- c(
    "ssm",
    class(out)
  )
  return(
    out
  )
}
