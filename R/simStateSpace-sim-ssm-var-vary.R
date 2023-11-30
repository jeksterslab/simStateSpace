#' Simulate Data from a Vector Autoregressive Model
#' using a State Space Model Parameterization
#' for n > 1 Individuals (Individual-Varying Parameters)
#'
#' This function simulates data from a vector autoregressive model
#' using a state space model parameterization
#' for `n > 1` individuals.
#' In this model,
#' the parameters can vary across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0_sqrt`,
#'   `alpha`,
#'   `beta`,
#'   `psi_sqrt`, or
#'   `gamma_eta`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMIVary
#' @inherit SimSSMFixed return
#' @inherit SimSSM references
#'
#' @examples
#' # prepare parameters
#' # In this example, beta varies across individuals
#' set.seed(42)
#' k <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- list(null_vec)
#' sigma0_sqrt <- list(iden_sqrt)
#' alpha <- list(null_vec)
#' beta <- list(
#'   diag(x = 0.1, nrow = k),
#'   diag(x = 0.2, nrow = k),
#'   diag(x = 0.3, nrow = k),
#'   diag(x = 0.4, nrow = k),
#'   diag(x = 0.5, nrow = k)
#' )
#' psi_sqrt <- list(iden_sqrt)
#' time <- 50
#' burn_in <- 0
#' gamma_eta <- list(0.10 * diag(k))
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
#' ssm <- SimSSMVARIVary(
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
#' ssm <- SimSSMVARIVary(
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
SimSSMVARIVary <- function(n,
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
      .SimSSM0VARIVary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_sqrt = rep(x = psi_sqrt, length.out = n),
        time = time,
        burn_in = burn_in
      )
    )
  } else {
    return(
      .SimSSM1VARIVary(
        n = n,
        mu0 = rep(x = mu0, length.out = n),
        sigma0_sqrt = rep(x = sigma0_sqrt, length.out = n),
        alpha = rep(x = alpha, length.out = n),
        beta = rep(x = beta, length.out = n),
        psi_sqrt = rep(x = psi_sqrt, length.out = n),
        gamma_eta = rep(x = gamma_eta, length.out = n),
        x = x,
        time = time,
        burn_in = burn_in
      )
    )
  }
}
