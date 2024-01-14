#' Simulate Data from an Ornstein–Uhlenbeck Model
#' using a State Space Model Parameterization
#' for n > 1 Individuals (Fixed Parameters)
#'
#' This function simulates data from an Ornstein–Uhlenbeck model
#' using a state space model parameterization
#' for `n > 1` individuals.
#' In this model,
#' the parameters are invariant across individuals.
#' See details for more information.
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
#'     \boldsymbol{\varepsilon}_{i, t}
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
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t}
#'   }
#'   where
#'   \eqn{\boldsymbol{\mu}}
#'   is the long-term mean or equilibrium level,
#'   \eqn{\boldsymbol{\Phi}}
#'   is the rate of mean reversion,
#'   determining how quickly the variable returns to its mean,
#'   \eqn{\boldsymbol{\Sigma}}
#'   is the matrix of volatility
#'   or randomness in the process, and
#'   \eqn{\mathrm{d}\boldsymbol{W}}
#'   is a Wiener process or Brownian motion,
#'   which represents random fluctuations.
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
#'     \boldsymbol{\varepsilon}_{i, t}
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
#'     \mathrm{d} \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t}
#'   }
#'   where
#'   \eqn{\mathbf{x}_{i, t}} is a vector of covariates
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}} is the coefficient matrix
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
#'     \boldsymbol{\Gamma}_{\mathbf{y}}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\varepsilon}_{i, t}
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
#'   \eqn{\boldsymbol{\Gamma}_{\mathbf{y}}} is the coefficient matrix
#'   linking the covariates to the observed variables.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \mathrm{d} \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\Phi}
#'     \left(
#'     \boldsymbol{\mu}
#'     -
#'     \boldsymbol{\eta}_{i, t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{i, t} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMOU
#' @inheritParams SimSSMFixed
#' @inherit SimSSMFixed return
#' @inherit SimSSMOU references
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' p <- k <- 2
#' iden <- diag(p)
#' n <- 5
#' mu0 <- c(-3.0, 1.5)
#' sigma0 <- iden
#' mu <- c(5.76, 5.18)
#' phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
#' sigma <- matrix(
#'   data = c(2.79, 0.06, 0.06, 3.27),
#'   nrow = p
#' )
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- diag(x = 0.50, nrow = k)
#' delta_t <- 0.10
#' time <- 50
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
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
#' ssm <- SimSSMOUFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   mu = mu,
#'   phi = phi,
#'   sigma = sigma,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   type = 0,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMOUFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   mu = mu,
#'   phi = phi,
#'   sigma = sigma,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMOUFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   mu = mu,
#'   phi = phi,
#'   sigma = sigma,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   delta_t = delta_t,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ou
#' @export
SimSSMOUFixed <- function(n,
                          mu0,
                          sigma0,
                          mu,
                          phi,
                          sigma,
                          nu,
                          lambda,
                          theta,
                          gamma_y = NULL,
                          gamma_eta = NULL,
                          x = NULL,
                          type = 0,
                          delta_t,
                          time,
                          burn_in = 0) {
  sigma0_l <- t(chol(sigma0))
  sigma_l <- t(chol(sigma))
  theta_l <- t(chol(theta))
  data <- switch(
    EXPR = as.character(type),
    "0" = {
      .SimSSM0OUFixed(
        n = n,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        mu = mu,
        phi = phi,
        sigma_l = sigma_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        delta_t = delta_t,
        time = time,
        burn_in = burn_in
      )
    },
    "1" = {
      .SimSSM1OUFixed(
        n = n,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        mu = mu,
        phi = phi,
        sigma_l = sigma_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        gamma_eta = gamma_eta,
        x = x,
        delta_t = delta_t,
        time = time,
        burn_in = burn_in
      )
    },
    "2" = {
      .SimSSM2OUFixed(
        n = n,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        mu = mu,
        phi = phi,
        sigma_l = sigma_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        gamma_y = gamma_y,
        gamma_eta = gamma_eta,
        x = x,
        delta_t = delta_t,
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
      mu = mu,
      phi = phi,
      sigma = sigma,
      nu = nu,
      lambda = lambda,
      theta = theta,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = type,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in,
      sigma0_l = sigma0_l,
      sigma_l = sigma_l,
      theta_l = theta_l
    ),
    model = list(
      model = "ou",
      n1 = FALSE,
      covariates = covariates,
      fixed = TRUE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMOUFixed"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
