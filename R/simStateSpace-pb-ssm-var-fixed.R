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
                          alpha = NULL, beta, psi_l,
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
  p <- dim(psi_l)[1]
  k <- p
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
  beta_label <- params_latent <- params_inicov <- matrix(
    data = "",
    nrow = p,
    ncol = p
  )
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      beta_label[i, j] <- paste0(
        "beta_",
        i,
        j
      )
      params_inicov[i, j] <- paste0("sigma0_", i, j)
      params_latent[i, j] <- paste0("psi_", i, j)
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
        "(beta_",
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
    beta
  )
  names(startval) <- c(
    beta_label
  )
  if (is.null(alpha)) {
    alpha_value <- rep(
      x = 0,
      time = p
    )
  } else {
    alpha_value <- alpha
    formula <- lapply(
      X = seq_len(length(formula)),
      FUN = function(x) {
        paste0(
          formula[[x]],
          " + ",
          "alpha_",
          x
        )
      }
    )
    startval <- c(
      beta,
      alpha
    )
    names(startval) <- c(
      beta_label,
      paste0(
        "alpha_",
        seq_len(p)
      )
    )
  }
  formula <- lapply(
    X = formula,
    FUN = stats::as.formula
  )
  y_names <- paste0("y", seq_len(k))
  eta_names <- paste0("eta_", seq_len(p))
  sigma0 <- tcrossprod(sigma0_l)
  psi <- tcrossprod(psi_l)
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
    values.int = NULL,
    params.int = NULL
  )
  dynr_dynamics <- dynr::prep.formulaDynamics(
    formula = formula,
    startval = startval,
    isContinuousTime = TRUE
  )
  dynr_noise <- dynr::prep.noise(
    values.latent = psi,
    params.latent = params_latent,
    values.observed = matrix(
      data = 0,
      nrow = k,
      ncol = k
    ),
    params.observed = matrix(
      data = "fixed",
      nrow = k,
      ncol = k
    )
  )
  foo <- function(i) {
    dynr_data <- dynr::dynr.data(
      dataframe = as.data.frame(
        x = SimSSMVARFixed(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          alpha = alpha_value,
          beta = beta,
          psi_l = psi_l,
          type = type,
          x = x,
          gamma = gamma
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
  psi_vec <- .Vech(psi)
  names(psi_vec) <- .Vech(params_latent)
  est <- c(
    startval,
    mu0_vec,
    sigma0_vec,
    psi_vec
  )
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
  # nolint end
}
