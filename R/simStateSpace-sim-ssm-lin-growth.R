#' Simulate Data from a Linear Growth Curve Model
#'
#' This function simulates data
#' from a linear growth curve model
#' for `n > 1` individuals.
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
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta
#'     \right)
#'   }
#'   where \eqn{y_{i, t}}, \eqn{\eta_{0_{i, t}}},
#'   \eqn{\eta_{1_{i, t}}},
#'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   are random variables and
#'   \eqn{\theta} is a model parameter.
#'   \eqn{y_{i, t}} is a vector of observed random variables
#'   at time \eqn{t} and individual \eqn{i},
#'   \eqn{\eta_{0_{i, t}}}
#'   and
#'   \eqn{\eta_{1_{i, t}}} form a vector of latent random variables
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#'   is a vector of random measurement errors
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\theta} is the variance of
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
#'     \eta_{0_{i, t - 1}} \\
#'     \eta_{1_{i, t - 1}} \\
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
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta
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
#'     \eta_{0_{i, t - 1}} \\
#'     \eta_{1_{i, t - 1}} \\
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
#'     \boldsymbol{\varepsilon}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \boldsymbol{\varepsilon}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     0,
#'     \theta
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
#'     \eta_{0_{i, t - 1}} \\
#'     \eta_{1_{i, t - 1}} \\
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
#' @param sigma0 Numeric matrix.
#'   The covariance matrix
#'   of the intercept and the slope.
#' @param theta Numeric.
#'   The common measurement error variance.
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
#'
#' @inheritParams SimSSMFixed
#' @inherit SimSSMFixed return
#' @inherit SimSSM references
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' n <- 10
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
#' theta <- 0.6
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
#'   sigma0 = sigma0,
#'   theta = theta,
#'   type = 0,
#'   time = time
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   theta = theta,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   time = time
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinGrowth(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   theta = theta,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   time = time
#' )
#'
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim growth
#' @export
SimSSMLinGrowth <- function(n,
                            mu0,
                            sigma0,
                            theta,
                            gamma_y = NULL,
                            gamma_eta = NULL,
                            x = NULL,
                            type = 0,
                            time) {
  sigma0_l <- t(chol(sigma0))
  theta_l <- sqrt(theta)
  data <- switch(
    EXPR = as.character(type),
    "0" = {
      .SimSSM0LinGrowth(
        n = n,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        theta_l = theta_l,
        time = time
      )
    },
    "1" = {
      stopifnot(
        !is.null(x),
        !is.null(gamma_eta)
      )
      .SimSSM1LinGrowth(
        n = n,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        theta_l = theta_l,
        gamma_eta = gamma_eta,
        x = x,
        time = time
      )
    },
    "2" = {
      stopifnot(
        !is.null(x),
        !is.null(gamma_y),
        !is.null(gamma_eta)
      )
      .SimSSM2LinGrowth(
        n = n,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        theta_l = theta_l,
        gamma_y = gamma_y,
        gamma_eta = gamma_eta,
        x = x,
        time = time
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
      theta = theta,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = type,
      time = time,
      sigma0_l = sigma0_l,
      theta_l = theta_l
    ),
    model = list(
      model = "lingrowth",
      n1 = FALSE,
      covariates = covariates,
      fixed = TRUE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMLinGrowth"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
