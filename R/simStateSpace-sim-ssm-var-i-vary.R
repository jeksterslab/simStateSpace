#' Simulate Data from the Vector Autoregressive Model
#' (Individual-Varying Parameters)
#'
#' This function simulates data from the
#' vector autoregressive model
#' using a state space model parameterization.
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
#'   `gamma`, or
#'   `kappa`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param type Integer.
#'   State space model type.
#'   See Details in [SimSSMVARFixed()] for more information.
#' @inheritParams SimSSMIVary
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
#'
#' # Type 0
#' ssm <- SimSSMVARIVary(
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
#' ssm <- SimSSMVARIVary(
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
SimSSMVARIVary <- function(n, time,
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
    data <- .SimSSMLatIVary0(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      alpha = rep(x = alpha, length.out = n),
      beta = rep(x = beta, length.out = n),
      psi_l = rep(x = psi_l, length.out = n)
    )
  }
  if (type == 1) {
    stopifnot(
      !is.null(x),
      !is.null(gamma)
    )
    data <- .SimSSMLatIVary1(
      n = n,
      time = time,
      delta_t = 1.0,
      mu0 = rep(x = mu0, length.out = n),
      sigma0_l = rep(x = sigma0_l, length.out = n),
      alpha = rep(x = alpha, length.out = n),
      beta = rep(x = beta, length.out = n),
      psi_l = rep(x = psi_l, length.out = n),
      x = x, gamma = rep(x = gamma, length.out = n)
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
      fixed = FALSE,
      vary_i = TRUE
    ),
    data = data,
    fun = "SimSSMVARIVary"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
