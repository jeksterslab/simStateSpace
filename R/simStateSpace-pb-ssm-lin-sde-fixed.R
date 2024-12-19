#' Parametric Bootstrap for the
#' Linear Stochastic Differential Equation Model
#' using a State Space Model Parameterization
#' (Fixed Parameters)
#'
#' This function simulates data from
#' a linear stochastic differential equation model
#' using a state-space model parameterization
#' and fits the model using the [dynr] package.
#' The process is repeated `R` times.
#' It assumes that the parameters remain constant
#' across individuals and over time.
#' At the momennt, the function only supports
#' `type = 0`.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams SimSSMLinSDEFixed
#' @inheritParams dynr::dynr.cook
#' @inherit SimSSMLinSDEFixed references details
#' @param R Positive integer.
#'   Number of bootstrap samples.
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
#' @param max_eval Positive integer.
#'   Maximum evaluation.
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
                             n, time, delta_t = 0.1,
                             mu0, sigma0_l,
                             iota, phi, sigma_l,
                             nu, lambda, theta_l,
                             type = 0,
                             x = NULL, gamma = NULL, kappa = NULL,
                             mu0_fixed = FALSE,
                             sigma0_fixed = FALSE,
                             alpha_level = 0.05,
                             max_eval = 100000,
                             optimization_flag = TRUE,
                             hessian_flag = FALSE,
                             verbose = FALSE,
                             weight_flag = FALSE,
                             debug_flag = FALSE,
                             perturb_flag = FALSE,
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
  p <- dim(lambda)[2]
  k <- dim(lambda)[1]
  # nocov start
  if (interactive()) {
    message(
      paste0(
        "\n",
        "Bootstrapping is computationally intensive.",
        "\n",
        "Consider using the argument `ncores` ",
        "to take advantage of multiple CPU cores.",
        "\n"
      )
    )
  }
  # nocov end
  # the function is limited to diagonal lambda for now
  covariates <- x
  args <- list(
    R = R,
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
    max_eval = max_eval,
    optimization_flag = optimization_flag,
    hessian_flag = hessian_flag,
    verbose = verbose,
    weight_flag = weight_flag,
    debug_flag = debug_flag,
    perturb_flag = perturb_flag,
    ncores = ncores,
    seed = seed
  )
  dynr_initial <- .DynrInitial(
    mu0 = mu0,
    sigma0_l = sigma0_l,
    mu0_fixed = mu0_fixed,
    sigma0_fixed = sigma0_fixed
  )
  mu0_values <- .Vec(
    dynr_initial$values.inistate[[1]]
  )
  mu0_labels <- .Vec(
    dynr_initial$params.inistate[[1]]
  )
  names(mu0_values) <- mu0_labels
  sigma0_values <- .Vech(
    dynr_initial$values.inicov[[1]]
  )
  sigma0_labels <- .Vech(
    dynr_initial$params.inicov[[1]]
  )
  names(sigma0_values) <- sigma0_labels
  dynr_measurement <- .DynrMeasurement(
    lambda = lambda,
    nu = nu
  )
  nu_values <- .Vec(
    dynr_measurement$values.int[[1]]
  )
  nu_labels <- .Vec(
    dynr_measurement$params.int[[1]]
  )
  names(nu_values) <- nu_labels
  dynr_noise <- .DynrNoise(
    process_l = sigma_l,
    theta_l = theta_l,
    continuous = FALSE
  )
  sigma_values <- .Vec(
    dynr_noise$values.latent[[1]]
  )
  sigma_labels <- .Vec(
    dynr_noise$params.latent[[1]]
  )
  names(sigma_values) <- sigma_labels
  theta_values <- .Vech(
    dynr_noise$values.observed[[1]]
  )
  theta_labels <- .Vech(
    dynr_noise$params.observed[[1]]
  )
  names(theta_values) <- theta_labels
  dynr_dynamics <- .DynrDynamics(
    dynamics = phi,
    intercept = iota,
    continuous = FALSE
  )
  dynamics_values <- dynr_dynamics$startval
  dynamics_labels <- names(dynamics_values)
  dynr_dynamics <- dynr_dynamics$dynamics
  est <- c(
    dynamics_values,
    sigma_values,
    nu_values,
    theta_values,
    mu0_values,
    sigma0_values
  )
  foo <- function(i,
                  n,
                  time,
                  delta_t,
                  mu0,
                  sigma0_l,
                  iota,
                  phi,
                  sigma_l,
                  nu,
                  lambda,
                  theta_l,
                  type,
                  covariates,
                  gamma,
                  kappa,
                  dynr_initial,
                  dynr_measurement,
                  dynr_noise,
                  dynr_dynamics,
                  max_eval,
                  optimization_flag,
                  hessian_flag,
                  verbose,
                  weight_flag,
                  debug_flag,
                  perturb_flag) {
    temp <- tempdir()
    outfile <- tempfile(
      pattern = "dynr_",
      tmpdir = temp,
      fileext = ".c"
    )
    on.exit(
      unlink(temp)
    )
    dynr_model <- dynr::dynr.model(
      data = .DynrData(
        object = SimSSMLinSDEFixed(
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
          x = covariates,
          gamma = gamma,
          kappa = kappa
        )
      ),
      initial = dynr_initial,
      measurement = dynr_measurement,
      dynamics = dynr_dynamics,
      noise = dynr_noise,
      outfile = outfile
    )
    dynr_model@options$maxeval <- max_eval
    if (verbose) {
      output <- dynr::dynr.cook(
        dynr_model,
        optimization_flag = optimization_flag,
        hessian_flag = hessian_flag,
        verbose = verbose,
        weight_flag = weight_flag,
        debug_flag = debug_flag,
        perturb_flag = perturb_flag
      )
    } else {
      utils::capture.output(
        output <- dynr::dynr.cook(
          dynr_model,
          optimization_flag = optimization_flag,
          hessian_flag = hessian_flag,
          verbose = verbose,
          weight_flag = weight_flag,
          debug_flag = debug_flag,
          perturb_flag = perturb_flag
        )
      )
    }
    return(
      output
    )
  }
  # nocov start
  par <- FALSE
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (ncores > 1) {
      par <- TRUE
    }
  }
  if (par) {
    cl <- parallel::makeCluster(ncores)
    on.exit(
      parallel::stopCluster(cl = cl)
    )
    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(
        cl = cl,
        iseed = seed
      )
    }
    fit <- parallel::parLapply(
      cl = cl,
      X = seq_len(R),
      fun = foo,
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
      dynr_initial = dynr_initial,
      dynr_measurement = dynr_measurement,
      dynr_noise = dynr_noise,
      dynr_dynamics = dynr_dynamics,
      max_eval = max_eval,
      optimization_flag = optimization_flag,
      hessian_flag = hessian_flag,
      verbose = verbose,
      weight_flag = weight_flag,
      debug_flag = debug_flag,
      perturb_flag = perturb_flag
    )
    thetahatstar <- parallel::parLapply(
      cl = cl,
      X = fit,
      fun = function(x) {
        x@transformed.parameters
      }
    )
  } else {
    # nocov end
    if (!is.null(seed)) {
      set.seed(seed)
    }
    fit <- lapply(
      X = seq_len(R),
      FUN = foo,
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
      dynr_initial = dynr_initial,
      dynr_measurement = dynr_measurement,
      dynr_noise = dynr_noise,
      dynr_dynamics = dynr_dynamics,
      max_eval = max_eval,
      optimization_flag = optimization_flag,
      hessian_flag = hessian_flag,
      verbose = verbose,
      weight_flag = weight_flag,
      debug_flag = debug_flag,
      perturb_flag = perturb_flag
    )
    thetahatstar <- lapply(
      X = fit,
      FUN = function(x) {
        x@transformed.parameters
      }
    )
  }
  out <- list(
    call = match.call(),
    args = args,
    fun = "PBSSMLinSDEFixed",
    thetahatstar = thetahatstar,
    vcov = stats::var(
      do.call(
        what = "rbind",
        args = thetahatstar
      )
    ),
    est = est[names(thetahatstar[[1]])],
    fit = fit
  )
  class(out) <- c(
    "statespacepb",
    class(out)
  )
  return(
    out
  )
}
