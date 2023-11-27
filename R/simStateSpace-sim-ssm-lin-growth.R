#' Simulate Data from a Linear Growth Curve Model
#'
#' This function simulates data
#' from a linear growth curve model.
#'
#' @details
#'   ## Type 0
#'
#'   The measurement model is given by
#'   \deqn{
#'     y_{i, t}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 0 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'       \eta_{0_{i, t}} \\
#'       \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta^{2}
#'     \right)
#'   }
#'   where \eqn{y_{i, t}}, \eqn{\eta_{0_{i, t}}},
#'   \eqn{\eta_{1_{i, t}}},
#'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   are random variables and
#'   and \eqn{\theta^{2}} is a model parameter.
#'   \eqn{y_{i, t}} is a vector of observed random variables
#'   at time \eqn{t} and individual \eqn{i},
#'   \eqn{\eta_{0_{i, t}}},
#'   \eqn{\eta_{1_{i, t}}} is forms a vector of latent random variables
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   is a vector of random measurement errors
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\theta^{2}} is the variance of
#'   \eqn{\boldsymbol{\varepsilon}}.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 1 \\
#'       0 & 1 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right) .
#'   }
#'
#'   The mean vector and covariance matrix of the intercept and slope
#'   are captured in the mean vector and covariance matrix
#'   of the initial condition given by
#'   \deqn{
#'     \boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}
#'     =
#'     \left(
#'     \begin{array}{c}
#'       \mu_{\eta_{0}} \\
#'       \mu_{\eta_{1}} \\
#'     \end{array}
#'     \right) \quad \mathrm{and,}
#'   }
#'
#'   \deqn{
#'     \boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       \sigma^{2}_{\eta_{0}} &
#'       \sigma_{\eta_{0}, \eta_{1}} \\
#'       \sigma_{\eta_{1}, \eta_{0}} &
#'       \sigma^{2}_{\eta_{1}} \\
#'     \end{array}
#'     \right) .
#'   }
#'
#'   ## Type 1
#'
#'   The measurement model is given by
#'   \deqn{
#'     y_{i, t}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 0 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'       \eta_{0_{i, t}} \\
#'       \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta^{2}
#'     \right) .
#'   }
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 1 \\
#'       0 & 1 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{i, t}
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
#'     y_{i, t}
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 0 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'       \eta_{0_{i, t}} \\
#'       \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
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
#'     0,
#'     \theta^{2}
#'     \right)
#'   }
#'   where
#'   \eqn{\boldsymbol{\Gamma}_{\mathbf{y}}} is the coefficient matrix
#'   linking the covariates to the observed variables.
#'
#'   The dynamic structure is given by
#'   \deqn{
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     =
#'     \left(
#'     \begin{array}{cc}
#'       1 & 1 \\
#'       0 & 1 \\
#'     \end{array}
#'     \right)
#'     \left(
#'     \begin{array}{c}
#'     \eta_{0_{i, t}} \\
#'     \eta_{1_{i, t}} \\
#'     \end{array}
#'     \right)
#'     +
#'     \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
#'     \mathbf{x}_{i, t} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu0 Numeric vector.
#'   A vector of length two.
#'   The first element is the mean of the intercept,
#'   and the second element is the mean of the slope.
#' @param sigma0_sqrt Numeric matrix.
#'   Cholesky decomposition of the covariance matrix
#'   of the intercept and the slope.
#' @param theta_sqrt Numeric.
#'   Square root of the common measurement error variance.
#' @param gamma_y Numeric matrix.
#'   Matrix relating the values of the covariate matrix
#'   at time `t` to `y` at time `t`
#'   (\eqn{\boldsymbol{\Gamma}_{\mathbf{y}}}).
#' @param gamma_eta Numeric matrix.
#'   Matrix relating the values of the covariate matrix
#'   at time `t` to the latent variables (intercept and slope) at time `t`
#'   (\eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}}).
#' @param x A list of length `n` of numeric matrices.
#'   Each element of the list
#'   is a matrix of observed covariates in `type = 1` or `type = 2`.
#'   The number of rows in each matrix should be equal to `time`.
#' @inheritParams SimSSMFixed
#' @inherit SimSSMFixed return
#' @inherit SimSSM references
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' n <- 5
#' mu0 <- c(0.615, 1.006)
#' sigma0 <- matrix(
#'   data = c(
#'     1.932,
#'     0.618,
#'     0.618,
#'     0.587
#'   ),
#'   nrow = 2
#' )
#' sigma0_sqrt <- chol(sigma0)
#' theta <- 0.6
#' theta_sqrt <- sqrt(theta)
#' time <- 10
#' gamma_y <- matrix(data = 0.10, nrow = 1, ncol = 2)
#' gamma_eta <- matrix(data = 0.10, nrow = 2, ncol = 2)
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = 2 * time),
#'         ncol = 2
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   theta_sqrt = theta_sqrt,
#'   type = 0,
#'   time = time
#' )
#'
#' str(ssm)
#'
#' # Type 1
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   theta_sqrt = theta_sqrt,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   time = time
#' )
#'
#' str(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   theta_sqrt = theta_sqrt,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   time = time
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim growth
#' @export
SimSSMLinGrowth <- function(n,
                            mu0,
                            sigma0_sqrt,
                            theta_sqrt,
                            gamma_y = NULL,
                            gamma_eta = NULL,
                            x = NULL,
                            type = 0,
                            time) {
  stopifnot(
    type %in% 0:2
  )
  if (type == 0) {
    return(
      .SimSSM0LinGrowth(
        n = n,
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        theta_sqrt = theta_sqrt,
        time = time
      )
    )
  }
  if (type == 1) {
    return(
      .SimSSM1LinGrowth(
        n = n,
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        theta_sqrt = theta_sqrt,
        gamma_eta = gamma_eta,
        x = x,
        time = time
      )
    )
  }
  if (type == 2) {
    return(
      .SimSSM2LinGrowth(
        n = n,
        mu0 = mu0,
        sigma0_sqrt = sigma0_sqrt,
        theta_sqrt = theta_sqrt,
        gamma_y = gamma_y,
        gamma_eta = gamma_eta,
        x = x,
        time = time
      )
    )
  }
}
