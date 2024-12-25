.PBSSMOUFixedData <- function(i,
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
                              kappa) {
  fn <- file.path(
    path,
    paste0(
      prefix,
      "_",
      "data",
      "_",
      i,
      ".Rds"
    )
  )
  run <- function(n,
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
                  kappa) {
    data <- SimSSMOUFixed(
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
      x = covariates,
      gamma = gamma,
      kappa = kappa
    )
    saveRDS(
      object = data,
      file = fn
    )
    invisible()
  }
  tryCatch(
    {
      if (file.exists(fn)) {
        data <- readRDS(
          file = fn
        )
      } else {
        run(
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
          kappa = kappa
        )
      }
    },
    error = function(e) {
      run(
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
        kappa = kappa
      )
    }
  )
  invisible()
}
