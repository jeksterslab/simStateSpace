#' Simulate Data from the Linear Stochastic Equation Model
#' using a State Space Model Parameterization (n = 1)
#'
#' This function simulates data from the linear stochastic equation model
#' using a state space model parameterization.
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
#'     \boldsymbol{\varepsilon}_{t},
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
#'     \mathrm{d} \boldsymbol{\eta}_{t}
#'     =
#'     \left(
#'     \mathbf{b}
#'     +
#'     \mathbf{A}
#'     \boldsymbol{\eta}_{t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \mathbf{Q}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{t}
#'   }
#'   where
#'   \eqn{\mathbf{b}}
#'   A term which is unobserved and constant over time,
#'   \eqn{\mathbf{A}}
#'   is the drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations,
#'   \eqn{\mathbf{Q}}
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
#'     \mathbf{y}_{t}
#'     =
#'     \boldsymbol{\nu}
#'     +
#'     \boldsymbol{\Lambda}
#'     \boldsymbol{\eta}_{t}
#'     +
#'     \boldsymbol{\varepsilon}_{t},
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
#'     \mathrm{d} \boldsymbol{\eta}_{t}
#'     =
#'     \left(
#'     \mathbf{b}
#'     +
#'     \mathbf{A}
#'     \boldsymbol{\eta}_{t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{t}
#'     +
#'     \mathbf{Q}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{t}
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
#'     \boldsymbol{\varepsilon}_{t},
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
#'     \mathrm{d} \boldsymbol{\eta}_{t}
#'     =
#'     \left(
#'     \mathbf{b}
#'     +
#'     \mathbf{A}
#'     \boldsymbol{\eta}_{t}
#'     \right)
#'     \mathrm{d}t
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{t}
#'     +
#'     \boldsymbol{\Sigma}^{\frac{1}{2}}
#'     \mathrm{d}
#'     \mathbf{W}_{t} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param b Numeric vector.
#'   A term which is unobserved and constant over time
#'   (\eqn{\mathbf{b}}).
#' @param a Numeric matrix.
#'   The drift matrix
#'   which represents the rate of change of the solution
#'   in the absence of any random fluctuations
#'   (\eqn{\mathbf{A}}).
#' @param q Numeric matrix.
#'   The matrix of volatility
#'   or randomness in the process
#'   (\eqn{\mathbf{Q}}).
#' @param delta_t Numeric.
#'   Time interval (\eqn{\delta_t}).
#' @inherit SimSSM return
#' @inheritParams SimSSM
#'
#' @references
#'   Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
#'   Continuous-time dynamic models: Connections to structural equation models
#'   and other discrete-time models.
#'   In R. H. Hoyle (Ed.), Handbook of structural equation modeling (2nd ed.).
#'   The Guilford Press.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' p <- k <- 2
#' iden <- diag(p)
#' mu0 <- c(-3.0, 1.5)
#' sigma0 <- iden
#' a <- matrix(
#'   data = c(
#'     -0.10,
#'     0.05,
#'     0.05,
#'     -0.10
#'   ),
#'   nrow = p
#' )
#' b <- c(0.317, 0.230)
#' q <- matrix(
#'   data = c(2.79, 0.06, 0.06, 3.27),
#'   nrow = p
#' )
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- diag(x = 0.50, nrow = k)
#' delta_t <- 0.10
#' time <- 1000
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
#' x <- matrix(
#'   data = rnorm(n = k * (time + burn_in)),
#'   ncol = k
#' )
#'
#'
#' # Type 0
#' ssm <- SimSSMLinSDE(
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   b = b,
#'   a = a,
#'   q = q,
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
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim linsde
#' @export
SimSSMLinSDE <- function(mu0,
                         sigma0,
                         b,
                         a,
                         q,
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
  q_l <- t(chol(q))
  theta_l <- t(chol(theta))
  data <- switch(
    EXPR = as.character(type),
    "0" = {
      .SimSSM0LinSDE(
        mu0 = mu0,
        sigma0_l = sigma0_l,
        b = b,
        a = a,
        q_l = q_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        delta_t = delta_t,
        time = time,
        burn_in = burn_in
      )
    },
    "1" = {
      stopifnot(
        !is.null(x),
        !is.null(gamma_eta)
      )
      .SimSSM1LinSDE(
        mu0 = mu0,
        sigma0_l = sigma0_l,
        b = b,
        a = a,
        q_l = q_l,
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
      stopifnot(
        !is.null(x),
        !is.null(gamma_y),
        !is.null(gamma_eta)
      )
      .SimSSM2LinSDE(
        mu0 = mu0,
        sigma0_l = sigma0_l,
        b = b,
        a = a,
        q_l = q_l,
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
      mu0 = mu0,
      sigma0 = sigma0,
      b = b,
      a = a,
      q = q,
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
      q_l = q_l,
      theta_l = theta_l
    ),
    model = list(
      model = "linsde",
      n1 = TRUE,
      covariates = covariates,
      fixed = FALSE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMLinSDE"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
