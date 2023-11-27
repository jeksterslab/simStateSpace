#' Simulate Data from the Vector Autoregressive Model
#' using a State Space Model Parameterization (n = 1)
#'
#' This function simulates data from the vector autoregressive model
#' using a state space model parameterization.
#' See details for more information.
#'
#' @details The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\eta}_{t} .
#'   }
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \boldsymbol{\eta}_{t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{t - 1}
#'     +
#'     \boldsymbol{\zeta}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\eta}_{t}},
#'   \eqn{\boldsymbol{\eta}_{t - 1}},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{t}}
#'   are random variables,
#'   and
#'   \eqn{\boldsymbol{\alpha}},
#'   \eqn{\boldsymbol{\beta}},
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   are model parameters.
#'   \eqn{\boldsymbol{\eta}_{t}}
#'   is a vector of latent variables
#'   at time \eqn{t},
#'   \eqn{\boldsymbol{\eta}_{t - 1}}
#'   is a vector of latent variables
#'   at time \eqn{t - 1},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{t}}
#'   is a vector of dynamic noise
#'   at time \eqn{t}.
#'   \eqn{\boldsymbol{\alpha}}
#'   is a vector of intercepts,
#'   \eqn{\boldsymbol{\beta}}
#'   is a matrix of autoregression
#'   and cross regression coefficients,
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   is the covariance matrix of
#'   \eqn{\boldsymbol{\zeta}_{t}}.
#'
#'   Note that when `gamma_eta` and `x` are not `NULL`,
#'   the dynamic structure is given by
#'   \deqn{
#'     \boldsymbol{\eta}_{t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{t - 1}
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{t}
#'     +
#'     \boldsymbol{\zeta}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\zeta}_{t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Psi}
#'     \right)
#'   }
#'   where
#'   \eqn{\mathbf{x}_{t}} is a vector of covariates at time \eqn{t},
#'   and \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}} is the coefficient matrix
#'   linking the covariates to the latent variables.
#'
#' @inheritParams SimSSM
#' @inherit SimSSM return
#' @inherit SimSSM references
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
#' mu0 <- null_vec
#' sigma0_sqrt <- iden_sqrt
#' alpha <- null_vec
#' beta <- diag(x = 0.5, nrow = k)
#' psi_sqrt <- iden_sqrt
#' time <- 50
#' burn_in <- 0
#' gamma_eta <- 0.10 * diag(k)
#' x <- matrix(
#'   data = rnorm(n = k * (time + burn_in)),
#'   ncol = k
#' )
#'
#' # No covariates
#' ssm <- SimSSMVAR(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' # With covariates
#' ssm <- SimSSMVAR(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim var
#' @export
SimSSMVAR <- function(mu0,
                      sigma0_sqrt,
                      alpha,
                      beta,
                      psi_sqrt,
                      gamma_eta = NULL,
                      x = NULL,
                      time = 0,
                      burn_in = 0) {
  if (is.null(gamma_eta) || is.null(x)) {
    return(
      .SimSSM0VAR(
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        alpha = alpha,
        beta = beta,
        psi_sqrt = psi_sqrt,
        time = time,
        burn_in = burn_in
      )
    )
  } else {
    return(
      .SimSSM1VAR(
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        alpha = alpha,
        beta = beta,
        psi_sqrt = psi_sqrt,
        gamma_eta = gamma_eta,
        x = x,
        time = time,
        burn_in = burn_in
      )
    )
  }
}
