# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Convert Parameters from the Ornstein–Uhlenbeck Model
#' to State Space Model Parameterization
#'
#' This function converts parameters from the Ornstein–Uhlenbeck model
#' to state space model parameterization.
#' See details for more information.
#'
#' @details The state space parameters
#'   as a function of the  Ornstein–Uhlenbeck model parameters
#'   are given by
#'   \deqn{
#'       \boldsymbol{\beta}
#'       =
#'       \exp{
#'         \left(
#'           - \boldsymbol{\Phi}
#'           \Delta_{t}
#'         \right)
#'       }
#'   }
#'
#'   \deqn{
#'       \boldsymbol{\alpha}
#'       =
#'       - \boldsymbol{\Phi}^{-1}
#'       \left(
#'         \boldsymbol{\beta} - \mathbf{I}_{p}
#'       \right)
#'   }
#'
#'   \deqn{
#'       \mathrm{vec}
#'       \left(
#'         \boldsymbol{\Psi}
#'       \right)
#'       =
#'       \left\{
#'         \left[
#'           \left(
#'             - \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
#'           \right)
#'           +
#'           \left(
#'             \mathbf{I}_{p} \otimes - \boldsymbol{\Phi}
#'           \right)
#'         \right]
#'         \left[
#'           \exp
#'           \left(
#'             \left[
#'               \left(
#'                 - \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
#'               \right)
#'               +
#'               \left(
#'                 \mathbf{I}_{p} \otimes - \boldsymbol{\Phi}
#'               \right)
#'             \right]
#'             \Delta_{t}
#'         \right)
#'         -
#'         \mathbf{I}_{p \times p}
#'       \right]
#'       \mathrm{vec}
#'       \left(
#'         \boldsymbol{\Sigma}
#'       \right)
#'     \right\}
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMOU
#'
#' @return Returns a list of state space parameters:
#'   - `alpha`: Numeric vector.
#'     Vector of intercepts for the dynamic model
#'     (\eqn{\boldsymbol{\alpha}}).
#'   - `beta`: Numeric matrix.
#'     Transition matrix relating the values of the latent variables
#'     at time `t - 1` to those at time `t`
#'     (\eqn{\boldsymbol{\beta}}).
#'   - `psi`: Numeric matrix.
#'     The process noise covariance matrix
#'     (\eqn{\boldsymbol{\Psi}}).
#'
#' @examples
#' p <- k <- 2
#' mu <- c(5.76, 5.18)
#' phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
#' sigma_sqrt <- chol(
#'   matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
#' )
#' delta_t <- 0.10
#'
#' OU2SSM(
#'   mu = mu,
#'   phi = phi,
#'   sigma_sqrt = sigma_sqrt,
#'   delta_t = delta_t
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim ou
#' @export
OU2SSM <- function(mu, phi, sigma_sqrt, delta_t) {
    .Call(`_simStateSpace_OU2SSM`, mu, phi, sigma_sqrt, delta_t)
}

.SimSSM0 <- function(mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0`, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, time, burn_in)
}

.SimSSM0Fixed <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0Fixed`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, time, burn_in)
}

.SimSSM0IVary <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0IVary`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, time, burn_in)
}

.SimSSM0LinGrowth <- function(n, mu0, sigma0_sqrt, theta_sqrt, time) {
    .Call(`_simStateSpace_SimSSM0LinGrowth`, n, mu0, sigma0_sqrt, theta_sqrt, time)
}

.SimSSM0LinGrowthIVary <- function(n, mu0, sigma0_sqrt, theta_sqrt, time) {
    .Call(`_simStateSpace_SimSSM0LinGrowthIVary`, n, mu0, sigma0_sqrt, theta_sqrt, time)
}

.SimSSM0OU <- function(mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0OU`, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, delta_t, time, burn_in)
}

.SimSSM0OUFixed <- function(n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0OUFixed`, n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, delta_t, time, burn_in)
}

.SimSSM0OUIVary <- function(n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0OUIVary`, n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, delta_t, time, burn_in)
}

.SimSSM0VAR <- function(mu0, sigma0_sqrt, alpha, beta, psi_sqrt, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0VAR`, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, time, burn_in)
}

.SimSSM0VARFixed <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0VARFixed`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, time, burn_in)
}

.SimSSM0VARIVary <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, time, burn_in) {
    .Call(`_simStateSpace_SimSSM0VARIVary`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, time, burn_in)
}

.SimSSM1 <- function(mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1`, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, time, burn_in)
}

.SimSSM1Fixed <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1Fixed`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, time, burn_in)
}

.SimSSM1IVary <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1IVary`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, time, burn_in)
}

.SimSSM1LinGrowth <- function(n, mu0, sigma0_sqrt, theta_sqrt, gamma_eta, x, time) {
    .Call(`_simStateSpace_SimSSM1LinGrowth`, n, mu0, sigma0_sqrt, theta_sqrt, gamma_eta, x, time)
}

.SimSSM1LinGrowthIVary <- function(n, mu0, sigma0_sqrt, theta_sqrt, gamma_eta, x, time) {
    .Call(`_simStateSpace_SimSSM1LinGrowthIVary`, n, mu0, sigma0_sqrt, theta_sqrt, gamma_eta, x, time)
}

.SimSSM1OU <- function(mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1OU`, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, delta_t, time, burn_in)
}

.SimSSM1OUFixed <- function(n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1OUFixed`, n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, delta_t, time, burn_in)
}

.SimSSM1OUIVary <- function(n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1OUIVary`, n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_eta, x, delta_t, time, burn_in)
}

.SimSSM1VAR <- function(mu0, sigma0_sqrt, alpha, beta, psi_sqrt, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1VAR`, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, gamma_eta, x, time, burn_in)
}

.SimSSM1VARFixed <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1VARFixed`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, gamma_eta, x, time, burn_in)
}

.SimSSM1VARIVary <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM1VARIVary`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, gamma_eta, x, time, burn_in)
}

.SimSSM2 <- function(mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM2`, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, time, burn_in)
}

.SimSSM2Fixed <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM2Fixed`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, time, burn_in)
}

.SimSSM2IVary <- function(n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, time, burn_in) {
    .Call(`_simStateSpace_SimSSM2IVary`, n, mu0, sigma0_sqrt, alpha, beta, psi_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, time, burn_in)
}

.SimSSM2LinGrowth <- function(n, mu0, sigma0_sqrt, theta_sqrt, gamma_y, gamma_eta, x, time) {
    .Call(`_simStateSpace_SimSSM2LinGrowth`, n, mu0, sigma0_sqrt, theta_sqrt, gamma_y, gamma_eta, x, time)
}

.SimSSM2LinGrowthIVary <- function(n, mu0, sigma0_sqrt, theta_sqrt, gamma_y, gamma_eta, x, time) {
    .Call(`_simStateSpace_SimSSM2LinGrowthIVary`, n, mu0, sigma0_sqrt, theta_sqrt, gamma_y, gamma_eta, x, time)
}

.SimSSM2OU <- function(mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM2OU`, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, delta_t, time, burn_in)
}

.SimSSM2OUFixed <- function(n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM2OUFixed`, n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, delta_t, time, burn_in)
}

.SimSSM2OUIVary <- function(n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, delta_t, time, burn_in) {
    .Call(`_simStateSpace_SimSSM2OUIVary`, n, mu0, sigma0_sqrt, mu, phi, sigma_sqrt, nu, lambda, theta_sqrt, gamma_y, gamma_eta, x, delta_t, time, burn_in)
}

