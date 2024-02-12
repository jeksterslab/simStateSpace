#' Simulate Data from the
#' Linear Growth Curve Model
#' (Individual-Varying Parameters)
#'
#' This function simulates data from the
#' linear growth curve model.
#' It assumes that the parameters can vary
#' across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0`,
#'   `mu`,
#'   `theta_l`,
#'   `gamma`, or
#'   `kappa`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @param type Integer.
#'   State space model type.
#'   See Details in [SimSSMLinGrowth()] for more information.
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu0 A list of numeric vectors.
#'   Each element of the list
#'   is a vector of length two.
#'   The first element is the mean of the intercept,
#'   and the second element is the mean of the slope.
#' @param sigma0_l A list of numeric matrices.
#'   Each element of the list is the
#'   Cholesky factorization (`t(chol(sigma0))`)
#'   of the covariance matrix
#'   of the intercept and the slope.
#' @param theta_l A list numeric values.
#'   Each element of the list
#'   is the square root of the common measurement error variance.
#' @inheritParams SimSSMIVary
#'
#' @inherit SimSSMFixed references return
#'
#' @examples
#' # prepare parameters
#' # In this example, the mean vector of the intercept and slope vary.
#' # Specifically,
#' # there are two sets of values representing two latent classes.
#' set.seed(42)
#' ## number of individuals
#' n <- 10
#' ## time points
#' time <- 5
#' ## dynamic structure
#' p <- 2
#' mu0_1 <- c(0.615, 1.006) # lower starting point, higher growth
#' mu0_2 <- c(1.000, 0.500) # higher starting point, lower growth
#' mu0 <- list(mu0_1, mu0_2)
#' sigma0 <- matrix(
#'   data = c(
#'     1.932,
#'     0.618,
#'     0.618,
#'     0.587
#'   ),
#'   nrow = p
#' )
#' sigma0_l <- list(t(chol(sigma0)))
#' ## measurement model
#' k <- 1
#' theta <- 0.50
#' theta_l <- list(sqrt(theta))
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
#' ssm <- SimSSMLinGrowthIVary(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#'
#' # Type 1
#' ssm <- SimSSMLinGrowthIVary(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' plot(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinGrowthIVary(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
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
#' @keywords simStateSpace sim growth
#' @export
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim growth
#' @export
SimSSMLinGrowthIVary <- function(n, time,
                                 mu0, sigma0_l, theta_l,
                                 type = 0,
                                 x = NULL, gamma = NULL, kappa = NULL) {
  stopifnot(type %in% c(0, 1, 2))
  p <- 2
  k <- 1
  stopifnot(
    length(mu0[[1]]) == p,
    dim(sigma0_l[[1]]) == c(p, p)
  )
  theta_l <- lapply(
    X = theta_l,
    FUN = as.matrix
  )
  stopifnot(
    dim(theta_l[[1]]) == c(k, k)
  )
  covariates <- FALSE
  if (type > 0) {
    covariates <- TRUE
  }
  alpha <- list(
    rep(x = 0, times = p)
  )
  beta <- list(
    matrix(
      data = c(1, 0, 1, 1),
      nrow = p
    )
  )
  psi_l <- list(
    matrix(
      data = 0,
      nrow = p,
      ncol = p
    )
  )
  nu <- list(
    rep(x = 0, times = k)
  )
  lambda <- list(
    matrix(
      data = c(1, 0),
      nrow = k
    )
  )
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
      model = "lingrowth",
      covariates = covariates,
      fixed = FALSE,
      vary_i = TRUE
    ),
    data = data,
    fun = "SimSSMLinGrowthIVary"
  )
  class(out) <- c(
    "simstatespace",
    class(out)
  )
  return(
    out
  )
}
