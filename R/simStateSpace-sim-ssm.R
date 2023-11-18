#' Simulate Data from a State Space Model (n = 1)
#'
#' This function simulates data from a state space model.
#' See details for more information.
#'
#' @details
#'   ## Type 0
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \boldsymbol{\Theta}
#'     \right)
#'   }
#'   where \eqn{\mathbf{y}_{t}}, \eqn{\boldsymbol{\eta}_{t}},
#'   and \eqn{\boldsymbol{\varepsilon}_{t}}
#'   are random variables and \eqn{\boldsymbol{\nu}},
#'   \eqn{\boldsymbol{\Lambda}},
#'   and \eqn{\boldsymbol{\Theta}} are model parameters.
#'   \eqn{\mathbf{y}_{t}} is a vector of observed random variables
#'   at time \eqn{t},
#'   \eqn{\boldsymbol{\eta}_{t}} is a vector of latent random variables
#'   at time \eqn{t},
#'   and \eqn{\boldsymbol{\varepsilon}_{t}}
#'   is a vector of random measurement errors
#'   at time \eqn{t},
#'   while \eqn{\boldsymbol{\nu}} is a vector of intercept,
#'   \eqn{\boldsymbol{\Lambda}} is a matrix of factor loadings,
#'   and \eqn{\boldsymbol{\Theta}} is the covariance matrix of
#'   \eqn{\boldsymbol{\varepsilon}}.
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
#'   where \eqn{\boldsymbol{\eta}_{t}}, \eqn{\boldsymbol{\eta}_{t - 1}},
#'   and \eqn{\boldsymbol{\zeta}_{t}} are random variables
#'   and \eqn{\boldsymbol{\alpha}}, \eqn{\boldsymbol{\beta}},
#'   and \eqn{\boldsymbol{\Psi}} are model parameters.
#'   \eqn{\boldsymbol{\eta}_{t}} is a vector of latent variables
#'   at time \eqn{t}, \eqn{\boldsymbol{\eta}_{t - 1}}
#'   is a vector of latent variables at
#'   time \eqn{t - 1},
#'   and \eqn{\boldsymbol{\zeta}_{t}} is a vector of dynamic noise
#'   at time \eqn{t} while \eqn{\boldsymbol{\alpha}}
#'   is a vector of intercepts,
#'   \eqn{\boldsymbol{\beta}} is a matrix of autoregression
#'   and cross regression coefficients,
#'   and \eqn{\boldsymbol{\Psi}} is the covariance matrix of
#'   \eqn{\boldsymbol{\zeta}_{t}}.
#'
#'   ## Type 1
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{t}
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
#'   ## Type 2
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\Gamma}_{\mathbf{y}}
#'     \boldsymbol{x}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{t}
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
#'     \boldsymbol{\eta}_{t}
#'     =
#'     \boldsymbol{\alpha}
#'     +
#'     \boldsymbol{\beta}
#'     \boldsymbol{\eta}_{t - 1}
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \boldsymbol{x}_{t}
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
#'     \right) .
#'   }
#'
#' @param mu0 Numeric vector.
#'   Mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#' @param sigma0_sqrt Numeric matrix.
#'   Cholesky decomposition of the covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#' @param alpha Numeric vector.
#'   Vector of intercepts for the dynamic model
#'   (\eqn{\boldsymbol{\alpha}}).
#' @param beta Numeric matrix.
#'   Transition matrix relating the values of the latent variables
#'   at time `t - 1` to those at time `t`
#'   (\eqn{\boldsymbol{\beta}}).
#' @param psi_sqrt Numeric matrix.
#'   Cholesky decomposition of the process noise covariance matrix
#'   (\eqn{\boldsymbol{\Psi}}).
#' @param nu Numeric vector.
#'   Vector of intercepts for the measurement model
#'   (\eqn{\boldsymbol{\nu}}).
#' @param lambda Numeric matrix.
#'   Factor loading matrix linking the latent variables
#'   to the observed variables
#'   (\eqn{\boldsymbol{\Lambda}}).
#' @param theta_sqrt Numeric matrix.
#'   Cholesky decomposition of the measurement error covariance matrix
#'   (\eqn{\boldsymbol{\Theta}}).
#' @param gamma_y Numeric matrix.
#'   Matrix relating the values of the covariate matrix
#'   at time `t` to the observed variables at time `t`
#'   (\eqn{\boldsymbol{\Gamma}_{\mathbf{y}}}).
#' @param gamma_eta Numeric matrix.
#'   Matrix relating the values of the covariate matrix
#'   at time `t` to the latent variables at time `t`
#'   (\eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}}).
#' @param x Numeric matrix.
#'   The matrix of observed covariates in `type = 1` or `type = 2`.
#'   The number of rows should be equal to `time + burn_in`.
#' @param type Integer.
#'   State space model type.
#' @param time Positive integer.
#'   Number of time points to simulate.
#' @param burn_in Positive integer.
#'   Number of burn-in points to exclude before returning the results.
#'
#' @references
#'   Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
#'   Equivalence and differences between structural equation modeling
#'   and state-space modeling techniques.
#'   *Structural Equation Modeling: A Multidisciplinary Journal*,
#'   17(2), 303–332.
#'   \doi{10.1080/10705511003661553}
#'
#' @return Returns a list with the following elements:
#'   - `y`: A `t` by `k` matrix of values for the manifest variables.
#'   - `eta`: A `t` by `p` matrix of values for the latent variables.
#'   - `time`: A vector of discrete time points from 1 to `t`.
#'   - `n`: Number of individuals.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
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
#' x <- matrix(
#'   data = rnorm(n = k * (time + burn_in)),
#'   ncol = k
#' )
#'
#' # Type 0
#' ssm <- SimSSM(
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
#' ssm <- SimSSM(
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
#' ssm <- SimSSM(
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
SimSSM <- function(mu0,
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
                   burn_in) {
  stopifnot(
    type %in% 0:2
  )
  if (type == 0) {
    return(
      .SimSSM0(
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
  }
  if (type == 1) {
    return(
      .SimSSM1(
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
  }
  if (type == 2) {
    return(
      .SimSSM2(
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
}
