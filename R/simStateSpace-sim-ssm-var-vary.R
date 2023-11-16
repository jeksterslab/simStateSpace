#' Simulate Data from a Vector Autoregressive Model
#' using a State Space Model Parameterization
#' for n > 1 Individuals (Varying Parameters)
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
#'   `beta`, and
#'   `psi_sqrt`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSM0Vary
#' @inherit SimSSM0Fixed return
#' @inherit SimSSM0 references
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
#'
#' ssm <- SimSSMVARVary(
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
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim var
#' @export
SimSSMVARVary <- function(n,
                          mu0,
                          sigma0_sqrt,
                          alpha,
                          beta,
                          psi_sqrt,
                          time,
                          burn_in) {
  stopifnot(
    is.list(mu0),
    is.list(sigma0_sqrt),
    is.list(alpha),
    is.list(beta),
    is.list(psi_sqrt)
  )
  stopifnot(
    length(mu0) <= n,
    length(sigma0_sqrt) <= n,
    length(alpha) <= n,
    length(beta) <= n,
    length(psi_sqrt) <= n
  )
  foo <- function(i,
                  mu0,
                  sigma0_sqrt,
                  alpha,
                  beta,
                  psi_sqrt) {
    dat_i <- SimSSMVAR(
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      time = time,
      burn_in = burn_in
    )
    return(
      list(
        y = dat_i$y,
        eta = dat_i$eta,
        time = dat_i$time,
        id = matrix(data = i, ncol = 1, nrow = time)
      )
    )
  }
  return(
    mapply(
      FUN = foo,
      i = seq_len(n),
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      SIMPLIFY = FALSE
    )
  )
}
