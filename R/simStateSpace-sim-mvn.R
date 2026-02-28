#' Simulate Vectors
#' from the Multivariate Normal Distribution
#'
#' This function simulates random vectors
#' from the multivariate normal distribution.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of replications.
#' @param mu Numeric vector.
#'   Mean vector (\eqn{\boldsymbol{\nu}}).
#' @param sigma_l Numeric matrix.
#'   Cholesky factorization (`t(chol(sigma))`)
#'   of the variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma}}.
#' @return Returns a list of random vectors.
#'
#' @examples
#' n <- 10
#' mu <- c(0, 0, 0)
#' sigma_l <- t(chol(0.001 * diag(3)))
#' SimMVN(n = n, mu = mu, sigma_l = sigma_l)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SimMVN <- function(n, mu, sigma_l) {
  .SimMVN(
    n = n,
    mu = mu,
    sigma_l = sigma_l
  )
}

#' Simulate Intercept Vectors
#' in a Discrete-Time Vector Autoregressive Model
#' from the Multivariate Normal Distribution
#'
#' This function simulates random intercept vectors
#' in a discrete-time vector autoregressive model
#' from the multivariate normal distribution.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of replications.
#' @param nu Numeric vector.
#'   Intercept (\eqn{\boldsymbol{\nu}}).
#' @param vcov_nu_l Numeric matrix.
#'   Cholesky factorization (`t(chol(vcov_nu))`)
#'   of the sampling variance-covariance matrix of
#'   \eqn{\boldsymbol{\nu}}.
#' @return Returns a list of random intercept vectors.
#'
#' @examples
#' n <- 10
#' nu <- c(0, 0, 0)
#' vcov_nu_l <- t(chol(0.001 * diag(3)))
#' SimNuN(n = n, nu = nu, vcov_nu_l = vcov_nu_l)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SimNuN <- function(n, nu, vcov_nu_l) {
  .SimMVN(
    n = n,
    mu = nu,
    sigma_l = vcov_nu_l
  )
}

#' Simulate Intercept Vectors
#' in a Discrete-Time Vector Autoregressive Model
#' from the Multivariate Normal Distribution
#'
#' This function simulates random intercept vectors
#' in a discrete-time vector autoregressive model
#' from the multivariate normal distribution.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of replications.
#' @param alpha Numeric vector.
#'   Intercept (\eqn{\boldsymbol{\alpha}}).
#' @param vcov_alpha_l Numeric matrix.
#'   Cholesky factorization (`t(chol(vcov_alpha))`)
#'   of the sampling variance-covariance matrix of
#'   \eqn{\boldsymbol{\alpha}}.
#' @return Returns a list of random intercept vectors.
#'
#' @examples
#' n <- 10
#' alpha <- c(0, 0, 0)
#' vcov_alpha_l <- t(chol(0.001 * diag(3)))
#' SimAlphaN(n = n, alpha = alpha, vcov_alpha_l = vcov_alpha_l)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SimAlphaN <- function(n, alpha, vcov_alpha_l) {
  .SimMVN(
    n = n,
    mu = alpha,
    sigma_l = vcov_alpha_l
  )
}

#' Simulate Intercept Vectors
#' in a Continuous-Time Vector Autoregressive Model
#' from the Multivariate Normal Distribution
#'
#' This function simulates random intercept vectors
#' in a continuous-time vector autoregressive model
#' from the multivariate normal distribution.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of replications.
#' @param iota Numeric vector.
#'   Intercept (\eqn{\boldsymbol{\iota}}).
#' @param vcov_iota_l Numeric matrix.
#'   Cholesky factorization (`t(chol(vcov_iota))`)
#'   of the sampling variance-covariance matrix of
#'   \eqn{\boldsymbol{\iota}}.
#' @return Returns a list of random intercept vectors.
#'
#' @examples
#' n <- 10
#' iota <- c(0, 0, 0)
#' vcov_iota_l <- t(chol(0.001 * diag(3)))
#' SimIotaN(n = n, iota = iota, vcov_iota_l = vcov_iota_l)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SimIotaN <- function(n, iota, vcov_iota_l) {
  .SimMVN(
    n = n,
    mu = iota,
    sigma_l = vcov_iota_l
  )
}

#' Simulate Set Point Vectors
#' from the Multivariate Normal Distribution
#'
#' This function simulates random set point vectors
#' from the multivariate normal distribution.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   Number of replications.
#' @param mu Numeric vector.
#'   Set point (\eqn{\boldsymbol{\mu}}).
#' @param vcov_mu_l Numeric matrix.
#'   Cholesky factorization (`t(chol(vcov_mu))`)
#'   of the sampling variance-covariance matrix of
#'   \eqn{\boldsymbol{\mu}}.
#' @return Returns a list of random set point vectors.
#'
#' @examples
#' n <- 10
#' mu <- c(0, 0, 0)
#' vcov_mu_l <- t(chol(0.001 * diag(3)))
#' SimMuN(n = n, mu = mu, vcov_mu_l = vcov_mu_l)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace ssm
#' @export
SimMuN <- function(n, mu, vcov_mu_l) {
  .SimMVN(
    n = n,
    mu = mu,
    sigma_l = vcov_mu_l
  )
}
