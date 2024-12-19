#' Parametric Bootstrap for the
#' Vector Autoregressive Model
#' (Fixed Parameters)
#'
#' This function simulates data from
#' a vector autoregressive model
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
#' @inheritParams SimSSMVARFixed
#' @inheritParams PBSSMLinSDEFixed
#' @inheritParams dynr::dynr.cook
#' @inherit SimSSMVARFixed references details
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
#'     \item{fun}{Function used ("PBSSMVARFixed").}
#'   }
#'
#' @examples
#' \dontrun{
#' # prepare parameters
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- 0.001 * diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- 0.001 * diag(p)
#' psi_l <- t(chol(psi))
#'
#' pb <- PBSSMVARFixed(
#'   R = 10L,
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
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
#' @keywords simStateSpace boot var
#' @export
PBSSMVARFixed <- function(R,
                          n, time,
                          mu0, sigma0_l,
                          alpha, beta, psi_l,
                          type = 0,
                          x = NULL, gamma = NULL,
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
  p <- dim(psi_l)[1]
  k <- p
  covariates <- x
  args <- list(
    R = R,
    n = n,
    time = time,
    mu0 = mu0,
    sigma0_l = sigma0_l,
    alpha = alpha,
    beta = beta,
    psi_l = psi_l,
    type = type,
    x = x,
    gamma = gamma,
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
    lambda = diag(p),
    nu = rep(
      x = 0,
      times = p
    )
  )
  dynr_noise <- .DynrNoise(
    process_l = psi_l,
    theta_l = matrix(
      data = 0,
      nrow = k,
      ncol = k
    ),
    continuous = FALSE
  )
  psi_values <- .Vec(
    dynr_noise$values.latent[[1]]
  )
  psi_labels <- .Vec(
    dynr_noise$params.latent[[1]]
  )
  names(psi_values) <- psi_labels
  theta_values <- .Vech(
    dynr_noise$values.observed[[1]]
  )
  theta_labels <- .Vech(
    dynr_noise$params.observed[[1]]
  )
  names(theta_values) <- theta_labels
  dynr_dynamics <- .DynrDynamics(
    dynamics = beta,
    intercept = alpha,
    continuous = FALSE
  )
  dynamics_values <- dynr_dynamics$startval
  dynr_dynamics <- dynr_dynamics$dynamics
  est <- c(
    dynamics_values,
    psi_values,
    theta_values,
    mu0_values,
    sigma0_values
  )
  foo <- function(i,
                  n,
                  time,
                  mu0,
                  sigma0_l,
                  alpha,
                  beta,
                  psi_l,
                  type,
                  covariates,
                  gamma,
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
        object = SimSSMVARFixed(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          alpha = alpha,
          beta = beta,
          psi_l = psi_l,
          type = type,
          x = covariates,
          gamma = gamma
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
      mu0 = mu0,
      sigma0_l = sigma0_l,
      alpha = alpha,
      beta = beta,
      psi_l = psi_l,
      type = type,
      covariates = covariates,
      gamma = gamma,
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
      mu0 = mu0,
      sigma0_l = sigma0_l,
      alpha = alpha,
      beta = beta,
      psi_l = psi_l,
      type = type,
      covariates = covariates,
      gamma = gamma,
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
    fun = "PBSSMVARFixed",
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
