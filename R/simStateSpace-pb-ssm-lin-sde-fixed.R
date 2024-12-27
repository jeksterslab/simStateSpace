#' Parametric Bootstrap for the
#' Linear Stochastic Differential Equation Model
#' using a State Space Model Parameterization
#' (Fixed Parameters)
#'
#' This function simulates data from
#' a linear stochastic differential equation model
#' using a state-space model parameterization
#' and fits the model using the `dynr` package.
#' The process is repeated `R` times.
#' It assumes that the parameters remain constant
#' across individuals and over time.
#' At the moment, the function only supports
#' `type = 0`.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMLinSDEFixed
#' @inheritParams dynr::dynr.cook
#' @inherit SimSSMLinSDEFixed references details
#' @param R Positive integer.
#'   Number of bootstrap samples.
#' @param path Path to a directory
#'   to store bootstrap samples and estimates.
#' @param prefix Character string.
#'   Prefix used for the file names
#'   for the bootstrap samples and estimates.
#' @param mu0_fixed Logical.
#'   If `mu0_fixed = TRUE`,
#'   fix the initial mean vector
#'   to `mu0`.
#'   If `mu0_fixed = FALSE`,
#'   `mu0` is estimated.
#' @param sigma0_fixed Logical.
#'   If `sigma0_fixed = TRUE`,
#'   fix the initial covariance matrix
#'   to `tcrossprod(sigma0_l)`.
#'   If `sigma0_fixed = FALSE`,
#'   `sigma0` is estimated.
#' @param alpha_level Numeric vector.
#'   Significance level \eqn{\alpha}.
#' @param ncores Positive integer.
#'   Number of cores to use.
#'   If `ncores = NULL`,
#'   use a single core.
#'   Consider using multiple cores
#'   when number of bootstrap samples `R`
#'   is a large value.
#' @param seed Random seed.
#' @param xtol_rel Stopping criteria option
#'   for parameter optimization.
#'   See [dynr::dynr.model()] for more details.
#' @param stopval Stopping criteria option
#'   for parameter optimization.
#'   See [dynr::dynr.model()] for more details.
#' @param ftol_rel Stopping criteria option
#'   for parameter optimization.
#'   See [dynr::dynr.model()] for more details.
#' @param ftol_abs Stopping criteria option
#'   for parameter optimization.
#'   See [dynr::dynr.model()] for more details.
#' @param maxeval Stopping criteria option
#'   for parameter optimization.
#'   See [dynr::dynr.model()] for more details.
#' @param maxtime Stopping criteria option
#'   for parameter optimization.
#'   See [dynr::dynr.model()] for more details.
#'
#' @return Returns an object
#'   of class `statespacepb` which is a list with the following elements:
#'   \describe{
#'     \item{call}{Function call.}
#'     \item{args}{Function arguments.}
#'     \item{thetahatstar}{Sampling distribution of
#'       \eqn{\boldsymbol{\hat{\theta}}}.}
#'     \item{vcov}{Sampling variance-covariance matrix of
#'       \eqn{\boldsymbol{\hat{\theta}}}.}
#'     \item{est}{Vector of estimated
#'       \eqn{\boldsymbol{\hat{\theta}}}.}
#'     \item{fun}{Function used ("PBSSMLinSDEFixed").}
#'   }
#'
#' @examples
#' \dontrun{
#' # prepare parameters
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' delta_t <- 0.10
#' ## dynamic structure
#' p <- 2
#' mu0 <- c(-3.0, 1.5)
#' sigma0 <- 0.001 * diag(p)
#' sigma0_l <- t(chol(sigma0))
#' iota <- c(0.317, 0.230)
#' phi <- matrix(
#'   data = c(
#'     -0.10,
#'     0.05,
#'     0.05,
#'     -0.10
#'   ),
#'   nrow = p
#' )
#' sigma <- matrix(
#'   data = c(
#'     2.79,
#'     0.06,
#'     0.06,
#'     3.27
#'   ),
#'   nrow = p
#' )
#' sigma_l <- t(chol(sigma))
#' ## measurement model
#' k <- 2
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.001 * diag(k)
#' theta_l <- t(chol(theta))
#'
#' pb <- PBSSMLinSDEFixed(
#'   R = 10L,
#'   path = getwd(),
#'   prefix = "lse",
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   iota = iota,
#'   phi = phi,
#'   sigma_l = sigma_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0,
#'   ncores = parallel::detectCores() - 1,
#'   seed = 42
#' )
#' print(pb)
#' summary(pb)
#' vcov(pb)
#' coef(pb)
#' confint(pb)
#' }
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace boot linsde
#' @export
PBSSMLinSDEFixed <- function(R,
                             path,
                             prefix,
                             n, time, delta_t = 0.1,
                             mu0, sigma0_l,
                             iota, phi, sigma_l,
                             nu, lambda, theta_l,
                             type = 0,
                             x = NULL, gamma = NULL, kappa = NULL,
                             mu0_fixed = FALSE,
                             sigma0_fixed = FALSE,
                             alpha_level = 0.05,
                             optimization_flag = TRUE,
                             hessian_flag = FALSE,
                             verbose = FALSE,
                             weight_flag = FALSE,
                             debug_flag = FALSE,
                             perturb_flag = FALSE,
                             xtol_rel = 1e-7,
                             stopval = -9999,
                             ftol_rel = -1,
                             ftol_abs = -1,
                             maxeval = as.integer(-1),
                             maxtime = -1,
                             ncores = NULL,
                             seed = NULL) {
  R <- as.integer(R)
  stopifnot(R > 0)
  if (!type == 0) {
    stop(
      paste0(
        "The function currently supports",
        "`type = 0`.",
        "\n"
      )
    )
  }
  covariates <- x
  args <- list(
    R = R,
    path = path,
    prefix = prefix,
    n = n,
    time = time,
    delta_t = delta_t,
    mu0 = mu0,
    sigma0_l = sigma0_l,
    iota = iota,
    phi = phi,
    sigma_l = sigma_l,
    nu = nu,
    lambda = lambda,
    theta_l = theta_l,
    type = type,
    x = x,
    gamma = gamma,
    kappa = kappa,
    mu0_fixed = mu0_fixed,
    sigma0_fixed = sigma0_fixed,
    alpha_level = alpha_level,
    optimization_flag = optimization_flag,
    hessian_flag = hessian_flag,
    verbose = verbose,
    weight_flag = weight_flag,
    debug_flag = debug_flag,
    perturb_flag = perturb_flag,
    xtol_rel = xtol_rel,
    stopval = stopval,
    ftol_rel = ftol_rel,
    ftol_abs = ftol_abs,
    maxeval = maxeval,
    maxtime = maxtime,
    ncores = ncores,
    seed = seed
  )
  par <- FALSE
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (ncores > 1) {
      par <- TRUE
    }
  }
  if (par) {
    os_type <- Sys.info()["sysname"]
    if (os_type == "Darwin") {
      fork <- TRUE
    } else if (os_type == "Linux") {
      fork <- TRUE
    } else {
      fork <- FALSE
    }
    if (fork) {
      output <- .PBSSMLinSDEFixedFork(
        R = R,
        path = path,
        prefix = prefix,
        n = n,
        time = time,
        delta_t = delta_t,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        iota = iota,
        phi = phi,
        sigma_l = sigma_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        type = type,
        covariates = covariates,
        gamma = gamma,
        kappa = kappa,
        optimization_flag = optimization_flag,
        hessian_flag = hessian_flag,
        verbose = verbose,
        weight_flag = weight_flag,
        debug_flag = debug_flag,
        perturb_flag = perturb_flag,
        xtol_rel = xtol_rel,
        stopval = stopval,
        ftol_rel = ftol_rel,
        ftol_abs = ftol_abs,
        maxeval = maxeval,
        maxtime = maxtime,
        ncores = ncores,
        seed = seed
      )
    } else {
      output <- .PBSSMLinSDEFixedSocket(
        R = R,
        path = path,
        prefix = prefix,
        n = n,
        time = time,
        delta_t = delta_t,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        iota = iota,
        phi = phi,
        sigma_l = sigma_l,
        nu = nu,
        lambda = lambda,
        theta_l = theta_l,
        type = type,
        covariates = covariates,
        gamma = gamma,
        kappa = kappa,
        optimization_flag = optimization_flag,
        hessian_flag = hessian_flag,
        verbose = verbose,
        weight_flag = weight_flag,
        debug_flag = debug_flag,
        perturb_flag = perturb_flag,
        xtol_rel = xtol_rel,
        stopval = stopval,
        ftol_rel = ftol_rel,
        ftol_abs = ftol_abs,
        maxeval = maxeval,
        maxtime = maxtime,
        ncores = ncores,
        seed = seed
      )
    }
  } else {
    output <- .PBSSMLinSDEFixedSerial(
      R = R,
      path = path,
      prefix = prefix,
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      iota = iota,
      phi = phi,
      sigma_l = sigma_l,
      nu = nu,
      lambda = lambda,
      theta_l = theta_l,
      type = type,
      covariates = covariates,
      gamma = gamma,
      kappa = kappa,
      optimization_flag = optimization_flag,
      hessian_flag = hessian_flag,
      verbose = verbose,
      weight_flag = weight_flag,
      debug_flag = debug_flag,
      perturb_flag = perturb_flag,
      xtol_rel = xtol_rel,
      stopval = stopval,
      ftol_rel = ftol_rel,
      ftol_abs = ftol_abs,
      maxeval = maxeval,
      maxtime = maxtime,
      seed = seed
    )
  }
  out <- list(
    call = match.call(),
    args = args,
    thetahatstar = output$thetahatstar,
    vcov = stats::var(
      do.call(
        what = "rbind",
        args = output$thetahatstar
      )
    ),
    est = output$prep$est[names(output$thetahatstar[[1]])],
    fun = "PBSSMLinSDEFixed"
  )
  class(out) <- c(
    "statespacepb",
    class(out)
  )
  return(
    out
  )
}
