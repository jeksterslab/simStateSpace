#' Simulate Data from a State Space Model (n = 1)
#'
#' This function simulates data from a state space model.
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
#'   where
#'   \eqn{\mathbf{y}_{t}},
#'   \eqn{\boldsymbol{\eta}_{t}},
#'   and
#'   \eqn{\boldsymbol{\varepsilon}_{t}}
#'   are random variables
#'   and
#'   \eqn{\boldsymbol{\nu}},
#'   \eqn{\boldsymbol{\Lambda}},
#'   and
#'   \eqn{\boldsymbol{\Theta}}
#'   are model parameters.
#'   \eqn{\mathbf{y}_{t}}
#'   is a vector of observed random variables,
#'   \eqn{\boldsymbol{\eta}_{t}}
#'   is a vector of latent random variables,
#'   and
#'   \eqn{\boldsymbol{\varepsilon}_{t}}
#'   is a vector of random measurement errors,
#'   at time \eqn{t}.
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
#'   \eqn{\mathbf{x}_{t}}
#'   is a vector of covariates
#'   at time \eqn{t},
#'   and
#'   \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}}
#'   is the coefficient matrix
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
#'     \mathbf{x}_{t}
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
#'   \eqn{\boldsymbol{\Gamma}_{\mathbf{y}}}
#'   is the coefficient matrix
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
#'     \right) .
#'   }
#'
#' @param mu0 Numeric vector.
#'   Mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#' @param sigma0 Numeric matrix.
#'   The covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#' @param alpha Numeric vector.
#'   Vector of intercepts for the dynamic model
#'   (\eqn{\boldsymbol{\alpha}}).
#' @param beta Numeric matrix.
#'   Transition matrix relating the values of the latent variables
#'   at time `t - 1` to those at time `t`
#'   (\eqn{\boldsymbol{\beta}}).
#' @param psi Numeric matrix.
#'   The process noise covariance matrix
#'   (\eqn{\boldsymbol{\Psi}}).
#' @param nu Numeric vector.
#'   Vector of intercepts for the measurement model
#'   (\eqn{\boldsymbol{\nu}}).
#' @param lambda Numeric matrix.
#'   Factor loading matrix linking the latent variables
#'   to the observed variables
#'   (\eqn{\boldsymbol{\Lambda}}).
#' @param theta Numeric matrix.
#'   The measurement error covariance matrix
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
#'   See Details for more information.
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
#' @return Returns an object of class `simstatespace`
#'   which is a list with the following elements:
#'   - `call`: Function call.
#'   - `args`: Function arguments.
#'   - `data`: Generated data which is a list of length `n`.
#'     `data` is a list with the following elements:
#'     * `id`: A vector of ones of length `t`.
#'     * `time`: A vector of time points of length `t`.
#'     * `y`: A `t` by `k` matrix of values for the manifest variables.
#'     * `eta`: A `t` by `p` matrix of values for the latent variables.
#'     * `x`: A `t` by `j` matrix of values for the covariates.
#'   - `fun`: Function used.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' null_vec <- rep(x = 0, times = k)
#' mu0 <- null_vec
#' sigma0 <- iden
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi <- iden
#' nu <- null_vec
#' lambda <- iden
#' theta <- diag(x = 0.50, nrow = k)
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
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSM(
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
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSM(
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
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ssm
#' @export
SimSSM <- function(mu0,
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
                   type = 0,
                   time,
                   burn_in = 0) {
  sigma0_l <- t(chol(sigma0))
  psi_l <- t(chol(psi))
  theta_l <- t(chol(theta))
  data <- switch(
    EXPR = as.character(type),
    "0" = {
      .SimSSM0(
        mu0 = mu0,
        sigma0_l = sigma0_l,
        alpha = alpha,
        beta = beta,
        psi_l = psi_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        time = time,
        burn_in = burn_in
      )
    },
    "1" = {
      .SimSSM1(
        mu0 = mu0,
        sigma0_l = sigma0_l,
        alpha = alpha,
        beta = beta,
        psi_l = psi_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        gamma_eta = gamma_eta,
        x = x,
        time = time,
        burn_in = burn_in
      )
    },
    "2" = {
      .SimSSM2(
        mu0 = mu0,
        sigma0_l = sigma0_l,
        alpha = alpha,
        beta = beta,
        psi_l = psi_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        gamma_y = gamma_y,
        gamma_eta = gamma_eta,
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
      n1 = TRUE,
      covariates = covariates,
      fixed = FALSE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSM"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
