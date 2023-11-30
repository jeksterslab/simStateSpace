#' Simulate Data using a State Space Model Parameterization
#' for n > 1 Individuals (Fixed Parameters)
#'
#' This function simulates data
#' using a state space model parameterization
#' for `n > 1` individuals.
#' In this model,
#' the parameters are invariant across individuals.
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
#'     \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{i, t - 1}
#'     +
#'     \boldsymbol{\zeta}_{i, t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\eta}_{i, t}},
#'   \eqn{\boldsymbol{\eta}_{i, t - 1}},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{i, t}}
#'   are random variables,
#'   and
#'   \eqn{\boldsymbol{\alpha}},
#'   \eqn{\boldsymbol{\beta}},
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   are model parameters.
#'   \eqn{\boldsymbol{\eta}_{i, t}}
#'   is a vector of latent variables
#'   at time \eqn{t} and individual \eqn{i},
#'   \eqn{\boldsymbol{\eta}_{i, t - 1}}
#'   is a vector of latent variables
#'   at time \eqn{t - 1} and individual \eqn{i},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{i, t}}
#'   is a vector of dynamic noise
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\boldsymbol{\alpha}}
#'   is a vector of intercepts,
#'   \eqn{\boldsymbol{\beta}}
#'   is a matrix of autoregression
#'   and cross regression coefficients,
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   is the covariance matrix of
#'   \eqn{\boldsymbol{\zeta}_{i, t}}.
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
#'     \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{i, t - 1}
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\zeta}_{i, t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right)
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
#'     \boldsymbol{\eta}_{i, t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{i, t - 1}
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\zeta}_{i, t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right) .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of individuals.
#' @param x A list of length `n` of numeric matrices.
#'   Each element of the list
#'   is a matrix of observed covariates in `type = 1` or `type = 2`.
#'   The number of rows in each matrix should be equal to `time + burn_in`.
#' @inheritParams SimSSM
#' @inherit SimSSM references
#'
#' @return Returns a list of length `n`.
#'   Each element is a list with the following elements:
#'   - `y`: A `t` by `k` matrix of values for the manifest variables.
#'   - `eta`: A `t` by `p` matrix of values for the latent variables.
#'   - `x`: A `t` by `j` matrix of values for the covariates.
#'   - `time`: A vector of discrete time points from 1 to `t`.
#'   - `id`: A vector of ID numbers of length `t`.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0_sqrt <- iden_sqrt
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi_sqrt <- iden_sqrt
#' nu <- null_vec
#' lambda <- iden
#' theta_sqrt <- chol(diag(x = 0.50, nrow = k))
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
#'    )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 1
#' ssm <- SimSSMFixed(
#'   n = n,
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
#'   type = 1,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # Type 2
#' ssm <- SimSSMFixed(
#'   n = n,
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
#'   type = 2,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ssm
#' @export
SimSSMFixed <- function(n,
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
                        type = 0,
                        time,
                        burn_in = 0) {
  stopifnot(
    type %in% 0:2
  )
  switch(
    EXPR = as.character(type),
    "0" = {
      return(
        .SimSSM0Fixed(
          n = n,
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
      )
    },
    "1" = {
      return(
        .SimSSM1Fixed(
          n = n,
          mu0 = mu0,
          sigma0_sqrt = sigma0_sqrt,
          alpha = alpha,
          beta = beta,
          psi_sqrt = psi_sqrt,
          nu = nu,
          lambda = lambda,
          theta_sqrt = theta_sqrt,
          gamma_eta = gamma_eta,
          x = x,
          time = time,
          burn_in = burn_in
        )
      )
    },
    "2" = {
      return(
        .SimSSM2Fixed(
          n = n,
          mu0 = mu0,
          sigma0_sqrt = sigma0_sqrt,
          alpha = alpha,
          beta = beta,
          psi_sqrt = psi_sqrt,
          nu = nu,
          lambda = lambda,
          theta_sqrt = theta_sqrt,
          gamma_y = gamma_y,
          gamma_eta = gamma_eta,
          x = x,
          time = time,
          burn_in = burn_in
        )
      )
    }
  )
}
