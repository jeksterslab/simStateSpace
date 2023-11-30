#' Simulate Data from a Vector Autoregressive Model
#' using a State Space Model Parameterization
#' for n > 1 Individuals (Fixed Parameters)
#'
#' This function simulates data from a vector autoregressive model
#' using a state space model parameterization
#' for `n > 1` individuals.
#' In this model,
#' the parameters are invariant across individuals.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMFixed
#' @inherit SimSSMFixed return
#' @inherit SimSSM references
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0_sqrt <- iden_sqrt
#' alpha <- null_vec
#' beta <- diag(x = 0.5, nrow = k)
#' psi_sqrt <- iden_sqrt
#' time <- 50
#' burn_in <- 0
#' gamma_eta <- 0.10 * diag(k)
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'     )
#'   }
#' )
#'
#' # No covariates
#' ssm <- SimSSMVARFixed(
#'   n = n,
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
#' ssm <- SimSSMVARFixed(
#'   n = n,
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
SimSSMVARFixed <- function(n,
                           mu0,
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
      .SimSSM0VARFixed(
        n = n,
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
      .SimSSM1VARFixed(
        n = n,
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
