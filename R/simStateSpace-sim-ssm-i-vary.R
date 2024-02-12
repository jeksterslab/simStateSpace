#' Simulate Data from a State Space Model
#' (Individual-Varying Parameters)
#'
#' This function simulates data using a state space model.
#' It assumes that the parameters can vary
#' across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0_l`,
#'   `alpha`,
#'   `beta`,
#'   `psi_l`,
#'   `nu`,
#'   `lambda`,
#'   `theta_l`,
#'   `gamma`, or
#'   `kappa`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMFixed
#' @param type Integer.
#'   State space model type.
#'   See Details in [SimSSMFixed()] for more information.
#' @param mu0 List of numeric vectors.
#'   Each element of the list
#'   is the mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#' @param sigma0_l List of numeric matrices.
#'   Each element of the list
#'   is the Cholesky factorization (`t(chol(sigma0))`)
#'   of the covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#' @param alpha List of numeric vectors.
#'   Each element of the list
#'   is the vector of constant values for the dynamic model
#'   (\eqn{\boldsymbol{\alpha}}).
#' @param beta List of numeric matrices.
#'   Each element of the list
#'   is the transition matrix relating the values of the latent variables
#'   at the previous to the current time point
#'   (\eqn{\boldsymbol{\beta}}).
#' @param psi_l List of numeric matrices.
#'   Each element of the list
#'   is the Cholesky factorization (`t(chol(psi))`)
#'   of the covariance matrix
#'   of the process noise
#'   (\eqn{\boldsymbol{\Psi}}).
#' @param nu List of numeric vectors.
#'   Each element of the list
#'   is the vector of intercept values for the measurement model
#'   (\eqn{\boldsymbol{\nu}}).
#' @param lambda List of numeric matrices.
#'   Each element of the list
#'   is the factor loading matrix linking the latent variables
#'   to the observed variables
#'   (\eqn{\boldsymbol{\Lambda}}).
#' @param theta_l List of numeric matrices.
#'   Each element of the list
#'   is the Cholesky factorization (`t(chol(theta))`)
#'   of the covariance matrix
#'   of the measurement error
#'   (\eqn{\boldsymbol{\Theta}}).
#' @param x List.
#'   Each element of the list is a matrix of covariates
#'   for each individual `i` in `n`.
#'   The number of columns in each matrix
#'   should be equal to `time`.
#' @param gamma List of numeric matrices.
#'   Each element of the list
#'   is the matrix linking the covariates to the latent variables
#'   at current time point
#'   (\eqn{\boldsymbol{\Gamma}}).
#' @param kappa List of numeric matrices.
#'   Each element of the list
#'   is the matrix linking the covariates to the observed variables
#'   at current time point
#'   (\eqn{\boldsymbol{\kappa}}).
#'
#' @inherit SimSSMFixed references return
#'
#' @examples
#' # prepare parameters
#' # In this example, beta varies across individuals.
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- list(
#'   rep(x = 0, times = p)
#' )
#' sigma0 <- 0.001 * diag(p)
#' sigma0_l <- list(
#'   t(chol(sigma0))
#' )
#' alpha <- list(
#'   rep(x = 0, times = p)
#' )
#' beta <- list(
#'   0.1 * diag(p),
#'   0.2 * diag(p),
#'   0.3 * diag(p),
#'   0.4 * diag(p),
#'   0.5 * diag(p)
#' )
#' psi <- 0.001 * diag(p)
#' psi_l <- list(
#'   t(chol(psi))
#' )
#' ## measurement model
#' k <- 3
#' nu <- list(
#'   rep(x = 0, times = k)
#' )
#' lambda <- list(
#'   diag(k)
#' )
#' theta <- 0.001 * diag(k)
#' theta_l <- list(
#'   t(chol(theta))
#' )
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
#' gamma <- list(
#'   diag(x = 0.10, nrow = p, ncol = j)
#' )
#' kappa <- list(
#'   diag(x = 0.10, nrow = k, ncol = j)
#' )
#'
#' # Type 0
#' ssm <- SimSSMIVary(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMIVary(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMIVary(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma = gamma,
#'   kappa = kappa
#' )
#'
#' plot(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ssm
#' @export
SimSSMIVary <- function(n, time, delta_t = 1.0,
                        mu0, sigma0_l,
                        alpha, beta, psi_l,
                        nu, lambda, theta_l,
                        type = 0,
                        x = NULL, gamma = NULL, kappa = NULL) {
  stopifnot(type %in% c(0, 1, 2))
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  if (type == 0) {
    data <- .SimSSMIVary0(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      alpha = rep(x = alpha, length.out = n),
      beta = rep(x = beta, length.out = n),
      psi_l = rep(x = psi_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n)
    )
  }
  if (type == 1) {
    stopifnot(
      !is.null(x),
      !is.null(gamma)
    )
    data <- .SimSSMIVary1(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      alpha = rep(x = alpha, length.out = n),
      beta = rep(x = beta, length.out = n),
      psi_l = rep(x = psi_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      x = rep(x = x, length.out = n),
      gamma = rep(x = gamma, length.out = n)
    )
  }
  if (type == 2) {
    stopifnot(
      !is.null(x),
      !is.null(gamma),
      !is.null(kappa)
    )
    data <- .SimSSMIVary2(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      alpha = rep(x = alpha, length.out = n),
      beta = rep(x = beta, length.out = n),
      psi_l = rep(x = psi_l, length.out = n),
      nu = rep(x = nu, length.out = n),
      lambda = rep(x = lambda, length.out = n),
      theta_l = rep(x = theta_l, length.out = n),
      x = rep(x = x, length.out = n),
      gamma = rep(x = gamma, length.out = n),
      kappa = rep(x = kappa, length.out = n)
    )
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n, time = time,
      mu0 = mu0, sigma0_l = sigma0_l,
      alpha = alpha, beta = beta, psi_l = psi_l,
      nu = nu, lambda = lambda, theta_l = theta_l,
      type = type,
      x = x, gamma = gamma, kappa = kappa
    ),
    model = list(
      model = "ssm",
      covariates = covariates,
      fixed = FALSE,
      vary_i = TRUE
    ),
    data = data,
    fun = "SimSSMIVary"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
