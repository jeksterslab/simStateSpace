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
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0 <- iden
#' alpha <- null_vec
#' beta <- diag(x = 0.5, nrow = k)
#' psi <- iden
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
#' SimSSMVARFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' # With covariates
#' SimSSMVARFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim var
#' @export
SimSSMVARFixed <- function(n,
                           mu0,
                           sigma0,
                           alpha,
                           beta,
                           psi,
                           gamma_eta = NULL,
                           x = NULL,
                           time = 0,
                           burn_in = 0) {
  sigma0_l <- t(chol(sigma0))
  psi_l <- t(chol(psi))
  if (is.null(gamma_eta) || is.null(x)) {
    data <- .SimSSM0VARFixed(
      n = n,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      alpha = alpha,
      beta = beta,
      psi_l = psi_l,
      time = time,
      burn_in = burn_in
    )
    covariates <- FALSE
  } else {
    data <- .SimSSM1VARFixed(
      n = n,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      alpha = alpha,
      beta = beta,
      psi_l = psi_l,
      gamma_eta = gamma_eta,
      x = x,
      time = time,
      burn_in = burn_in
    )
    covariates <- TRUE
  }
  out <- list(
    call = match.call(),
    args = list(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      gamma_eta = gamma_eta,
      x = x,
      time = time,
      burn_in = burn_in,
      sigma0_l = sigma0_l,
      psi_l = psi_l
    ),
    model = list(
      model = "var",
      n1 = FALSE,
      covariates = covariates,
      fixed = TRUE,
      vary_i = FALSE
    ),
    data = data,
    fun = "SimSSMVARFixed"
  )
  class(out) <- c(
    "ssm",
    class(out)
  )
  return(
    out
  )
}
