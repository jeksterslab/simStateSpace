#' Parametric Bootstrap for the
#' Ornstein–Uhlenbeck Model
#' using a State Space Model Parameterization
#' (Fixed Parameters)
#'
#' This function simulates data from
#' a Ornstein–Uhlenbeck (OU) model
#' using a state-space model parameterization
#' and fits the model using the [dynr] package.
#' The process is repeated `R` times.
#' It assumes that the parameters remain constant
#' across individuals and over time.
#' At the moment, the function only supports
#' `type = 0`.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMOUFixed
#' @inheritParams PBSSMLinSDEFixed
#' @inheritParams dynr::dynr.cook
#' @inherit SimSSMOUFixed references details
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
#'     \item{fun}{Function used ("PBSSMOUFixed").}
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
#' mu <- c(5.76, 5.18)
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
#' pb <- PBSSMOUFixed(
#'   R = 10L,
#'   path = getwd(),
#'   prefix = "ou",
#'   n = n,
#'   time = time,
#'   delta_t = delta_t,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   mu = mu,
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
#' @keywords simStateSpace boot ou
#' @export
PBSSMOUFixed <- function(R,
                         path,
                         prefix,
                         n, time, delta_t = 0.1,
                         mu0, sigma0_l,
                         mu, phi, sigma_l,
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
    mu = mu,
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
      output <- .PBSSMOUFixedFork(
        R = R,
        path = path,
        prefix = prefix,
        n = n,
        time = time,
        delta_t = delta_t,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        mu = mu,
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
      output <- .PBSSMOUFixedSocket(
        R = R,
        path = path,
        prefix = prefix,
        n = n,
        time = time,
        delta_t = delta_t,
        mu0 = mu0,
        sigma0_l = sigma0_l,
        mu = mu,
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
    output <- .PBSSMOUFixedSerial(
      R = R,
      path = path,
      prefix = prefix,
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      mu = mu,
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
    fun = "PBSSMOUFixed"
  )
  class(out) <- c(
    "statespacepb",
    class(out)
  )
  return(
    out
  )
}
