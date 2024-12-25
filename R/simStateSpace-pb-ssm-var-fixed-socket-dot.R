.PBSSMVARFixedSocket <- function(R,
                                 path,
                                 prefix,
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
  prep <- .PBSSMVARFixedPrepDynr(
    mu0 = mu0,
    sigma0_l = sigma0_l,
    alpha = alpha,
    beta = beta,
    psi_l = psi_l,
    mu0_fixed = mu0_fixed,
    sigma0_fixed = sigma0_fixed
  )
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
  parallel::clusterEvalQ(
    cl = cl,
    expr = library(simStateSpace)
  )
  parallel::clusterEvalQ(
    cl = cl,
    expr = library(dynr)
  )
  if (interactive()) {
    message(
      "\nGenerating data...\n"
    )
  }
  parallel::parLapply(
    cl = cl,
    X = seq_len(R),
    fun = .PBSSMVARFixedData,
    path = path,
    prefix = prefix,
    n = n,
    time = time,
    mu0 = mu0,
    sigma0_l = sigma0_l,
    alpha = alpha,
    beta = beta,
    psi_l = psi_l,
    type = type,
    covariates = covariates,
    gamma = gamma
  )
  if (interactive()) {
    message(
      "Model fitting...\n"
    )
  }
  parallel::parLapply(
    cl = cl,
    X = seq_len(R),
    fun = .PBFitDynr,
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
    maxtime = maxtime
  )
  thetahatstar <- parallel::parLapply(
    cl = cl,
    X = seq_len(R),
    FUN = .PBCoefDynr,
    path = path,
    prefix = prefix
  )
  return(
    list(
      prep = prep,
      thetahatstar = thetahatstar
    )
  )
}
