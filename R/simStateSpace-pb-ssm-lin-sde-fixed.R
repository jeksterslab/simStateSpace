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
                             iota = NULL, phi, sigma_l,
                             nu = NULL, theta_l,
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
  # nolint start
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
  p <- dim(sigma_l)[1]
  k <- dim(theta_l)[1]
  if (k != p) {
    stop(
      paste0(
        "The function currently supports",
        "single-indicator latent variable models.",
        "\n"
      )
    )
  }
  # the function is limited to diagonal lambda for now
  lambda <- diag(p)
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
  phi_label <- params_latent <- params_inicov <- matrix(
    data = "",
    nrow = p,
    ncol = p
  )
  params_observed <- matrix(
    data = "fixed",
    nrow = k,
    ncol = k
  )
  diag(params_observed) <- paste0(
    "theta_",
    seq_len(k),
    seq_len(k)
  )
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      phi_label[i, j] <- paste0(
        "phi_",
        i,
        j
      )
      params_inicov[i, j] <- paste0("sigma0_", i, j)
      params_latent[i, j] <- paste0("sigma_", i, j)
    }
  }
  params_inicov[
    upper.tri(params_inicov)
  ] <- t(params_inicov)[
    upper.tri(params_inicov)
  ]
  params_latent[
    upper.tri(params_latent)
  ] <- t(params_latent)[
    upper.tri(params_latent)
  ]
  params_inistate <- paste0(
    "mu0_",
    seq_len(p)
  )
  if (mu0_fixed) {
    params_inistate <- rep(
      x = "fixed",
      times = p
    )
  }
  if (sigma0_fixed) {
    params_inicov <- matrix(
      data = "fixed",
      nrow = p,
      ncol = p
    )
  }
  formula <- lapply(
    X = seq_len(p),
    FUN = function(i) {
      terms <- paste0(
        "(phi_",
        i,
        seq_len(p),
        " * eta_",
        seq_len(p),
        ")",
        collapse = " + "
      )
      return(
        paste0(
          "eta_",
          i,
          " ~ ",
          terms
        )
      )
    }
  )
  startval <- c(
    phi
  )
  names(startval) <- c(
    phi_label
  )
  if (is.null(iota)) {
    iota_value <- rep(
      x = 0,
      time = p
    )
  } else {
    iota_value <- iota
    formula <- lapply(
      X = seq_len(length(formula)),
      FUN = function(x) {
        paste0(
          formula[[x]],
          " + ",
          "iota_",
          x
        )
      }
    )
    startval <- c(
      phi,
      iota
    )
    names(startval) <- c(
      phi_label,
      paste0(
        "iota_",
        seq_len(p)
      )
    )
  }
  formula <- lapply(
    X = formula,
    FUN = stats::as.formula
  )
  if (is.null(nu)) {
    nu_value <- rep(
      x = 0,
      time = k
    )
    values_int <- NULL
    params_int <- NULL
  } else {
    nu_value <- nu
    values_int <- matrix(
      data = nu_value,
      ncol = 1
    )
    params_int <- matrix(
      data = paste0(
        "nu_",
        seq_len(k)
      ),
      ncol = 1
    )
  }
  y_names <- paste0("y", seq_len(k))
  eta_names <- paste0("eta_", seq_len(p))
  sigma0 <- tcrossprod(sigma0_l)
  sigma <- tcrossprod(sigma_l)
  theta <- tcrossprod(theta_l)
  dynr_initial <- dynr::prep.initial(
    values.inistate = mu0,
    params.inistate = params_inistate,
    values.inicov = sigma0,
    params.inicov = params_inicov
  )
  dynr_measurement <- dynr::prep.measurement(
    values.load = diag(p),
    params.load = matrix(
      data = "fixed",
      nrow = p,
      ncol = p
    ),
    state.names = eta_names,
    obs.names = y_names,
    values.int = values_int,
    params.int = params_int
  )
  dynr_dynamics <- dynr::prep.formulaDynamics(
    formula = formula,
    startval = startval,
    isContinuousTime = TRUE
  )
  dynr_noise <- dynr::prep.noise(
    values.latent = sigma,
    params.latent = params_latent,
    values.observed = theta,
    params.observed = params_observed
  )
  foo <- function(i) {
    dynr_data <- dynr::dynr.data(
      dataframe = as.data.frame(
        x = SimSSMLinSDEFixed(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota_value,
          phi = phi,
          sigma_l = sigma_l,
          nu = nu_value,
          lambda = lambda,
          theta_l = theta_l,
          type = type,
          x = x,
          gamma = gamma,
          kappa = kappa
        )
      ),
      id = "id",
      time = "time",
      observed = y_names
    )
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
      data = dynr_data,
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
    vars_to_export <- ls(environment())
    parallel::clusterExport(
      cl = cl,
      varlist = vars_to_export,
      envir = environment()
    )
    fit <- parallel::parLapply(
      cl = cl,
      X = seq_len(R),
      fun = foo
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
      FUN = foo
    )
    thetahatstar <- lapply(
      X = fit,
      FUN = function(x) {
        x@transformed.parameters
      }
    )
  }
  mu0_vec <- mu0
  names(mu0_vec) <- params_inistate
  sigma0_vec <- .Vech(sigma0)
  names(sigma0_vec) <- .Vech(params_inicov)
  sigma_vec <- .Vech(sigma)
  names(sigma_vec) <- .Vech(params_latent)
  theta_vec <- diag(theta)
  names(theta_vec) <- diag(params_observed)
  nu_vec <- nu_value
  names(nu_vec) <- paste0(
    "nu_",
    seq_len(k)
  )
  est <- c(
    startval,
    mu0_vec,
    sigma0_vec,
    sigma_vec,
    theta_vec,
    nu_vec
  )
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
  # nolint end
}
