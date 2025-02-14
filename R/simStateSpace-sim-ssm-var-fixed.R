#' Simulate Data from the Vector Autoregressive Model
#' (Fixed Parameters)
#'
#' This function simulates data from the
#' vector autoregressive model
#' using a state space model parameterization.
#' It assumes that the parameters remain constant
#' across individuals and over time.
#'
#' @details
#'   ## Type 0
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{i, t}
#'     =
#'     \boldsymbol{\eta}_{i, t}
#'   }
#'   where \eqn{\mathbf{y}_{i, t}}
#'   represents a vector of observed variables
#'   and \eqn{\boldsymbol{\eta}_{i, t}}
#'   a vector of latent variables
#'   for individual \eqn{i} and time \eqn{t}.
#'   Since the observed and latent variables are equal,
#'   we only generate data
#'   from the dynamic structure.
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
#'     \boldsymbol{\zeta}_{i, t},
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
#'   Here,
#'   \eqn{\boldsymbol{\eta}_{i, t}}
#'   is a vector of latent variables
#'   at time \eqn{t} and individual \eqn{i},
#'   \eqn{\boldsymbol{\eta}_{i, t - 1}}
#'   represents a vector of latent variables
#'   at time \eqn{t - 1} and individual \eqn{i},
#'   and
#'   \eqn{\boldsymbol{\zeta}_{i, t}}
#'   represents a vector of dynamic noise
#'   at time \eqn{t} and individual \eqn{i}.
#'   \eqn{\boldsymbol{\alpha}}
#'   denotes a vector of intercepts,
#'   \eqn{\boldsymbol{\beta}}
#'   a matrix of autoregression
#'   and cross regression coefficients,
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   the covariance matrix of
#'   \eqn{\boldsymbol{\zeta}_{i, t}}.
#'
#'   An alternative representation of the dynamic noise
#'   is given by
#'   \deqn{
#'     \boldsymbol{\zeta}_{i, t}
#'     =
#'     \boldsymbol{\Psi}^{\frac{1}{2}}
#'     \mathbf{z}_{i, t},
#'     \quad
#'     \mathrm{with}
#'     \quad
#'     \mathbf{z}_{i, t}
#'     \sim
#'     \mathcal{N}
#'     \left(
#'     \mathbf{0},
#'     \mathbf{I}
#'     \right)
#'   }
#'   where
#'   \eqn{
#'     \left( \boldsymbol{\Psi}^{\frac{1}{2}} \right)
#'     \left( \boldsymbol{\Psi}^{\frac{1}{2}} \right)^{\prime}
#'     =
#'     \boldsymbol{\Psi} .
#'   }
#'
#'   ## Type 1
#'
#'   The measurement model is given by
#'   \deqn{
#'     \mathbf{y}_{i, t}
#'     =
#'     \boldsymbol{\eta}_{i, t} .
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
#'     \boldsymbol{\Gamma}
#'     \mathbf{x}_{i, t}
#'     +
#'     \boldsymbol{\zeta}_{i, t},
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
#'   \eqn{\mathbf{x}_{i, t}} represents a vector of covariates
#'   at time \eqn{t} and individual \eqn{i},
#'   and \eqn{\boldsymbol{\Gamma}} the coefficient matrix
#'   linking the covariates to the latent variables.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMFixed
#'
#' @inherit SimSSMFixed references return
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- 0.001 * diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- 0.001 * diag(p)
#' psi_l <- t(chol(psi))
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     matrix(
#'       data = stats::rnorm(n = time * j),
#'       nrow = j,
#'       ncol = time
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMVARFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMVARFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim var
#' @export
SimSSMVARFixed <- function(n, time,
                           mu0, sigma0_l,
                           alpha, beta, psi_l,
                           type = 0,
                           x = NULL, gamma = NULL) {
  stopifnot(type %in% c(0, 1))
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  if (type == 0) {
    data <- .SimSSMLatFixed0(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l
    )
  }
  if (type == 1) {
    stopifnot(
      !is.null(x),
      !is.null(gamma)
    )
    data <- .SimSSMLatFixed1(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      x = x, gamma = gamma
    )
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n, time = time,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      type = type,
      x = x, gamma = gamma
    ),
    model = list(
      model = "var",
      covariates = covariates,
      fixed = TRUE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMVARFixed"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  out
}
