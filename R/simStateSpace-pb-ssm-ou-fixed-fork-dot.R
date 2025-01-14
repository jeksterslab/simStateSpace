.PBSSMOUFixedFork <- function(R,
                              path,
                              prefix,
                              n,
                              time,
                              delta_t,
                              mu0,
                              sigma0_l,
                              mu,
                              phi,
                              sigma_l,
                              nu,
                              lambda,
                              theta_l,
                              type,
                              covariates,
                              gamma,
                              kappa,
                              mu0_fixed,
                              sigma0_fixed,
                              optimization_flag,
                              hessian_flag,
                              verbose,
                              weight_flag,
                              debug_flag,
                              perturb_flag,
                              xtol_rel,
                              stopval,
                              ftol_rel,
                              ftol_abs,
                              maxeval,
                              maxtime,
                              ncores,
                              seed) {
  if (ncores > R) {
    ncores <- R
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  prep <- .PBSSMOUFixedPrepDynr(
    mu0 = mu0,
    sigma0_l = sigma0_l,
    mu = mu,
    phi = phi,
    sigma_l = sigma_l,
    nu = nu,
    lambda = lambda,
    theta_l = theta_l,
    mu0_fixed = mu0_fixed,
    sigma0_fixed = sigma0_fixed
  )
  if (interactive()) {
    message(
      "\nGenerating data...\n"
    )
  }
  parallel::mclapply(
    X = seq_len(R),
    FUN = .PBSSMOUFixedData,
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
    mc.cores = ncores
  )
  if (interactive()) {
    message(
      "Model fitting...\n"
    )
  }
  parallel::mclapply(
    X = seq_len(R),
    FUN = .PBFitDynr,
    path = path,
    prefix = prefix,
    dynr_initial = prep$dynr_initial,
    dynr_measurement = prep$dynr_measurement,
    dynr_noise = prep$dynr_noise,
    dynr_dynamics = prep$dynr_dynamics,
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
    mc.cores = ncores
  )
  thetahatstar <- parallel::mclapply(
    X = seq_len(R),
    FUN = .PBCoefDynr,
    path = path,
    prefix = prefix,
    mc.cores = ncores
  )
  return(
    list(
      prep = prep,
      thetahatstar = thetahatstar
    )
  )
}
