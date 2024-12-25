.PBSSMVARFixedData <- function(i,
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
                               gamma) {
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
                  mu0,
                  sigma0_l,
                  alpha,
                  beta,
                  psi_l,
                  type,
                  covariates,
                  gamma) {
    data <- SimSSMVARFixed(
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
          mu0 = mu0,
          sigma0_l = sigma0_l,
          alpha = alpha,
          beta = beta,
          psi_l = psi_l,
          type = type,
          covariates = covariates,
          gamma = gamma
        )
      }
    },
    error = function(e) {
      run(
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
    }
  )
  invisible()
}
