.DynrInitial <- function(mu0,
                         sigma0_l,
                         mu0_fixed,
                         sigma0_fixed) {
  sigma0 <- tcrossprod(sigma0_l)
  if (mu0_fixed) {
    params_inistate <- rep(
      x = "fixed",
      times = length(mu0)
    )
  } else {
    params_inistate <- paste0(
      "mu0_",
      seq_len(
        length(mu0)
      ),
      "_1"
    )
  }
  if (sigma0_fixed) {
    params_inicov <- matrix(
      data = "fixed",
      nrow = dim(sigma0)[1],
      ncol = dim(sigma0)[2]
    )
  } else {
    params_inicov <- .LabelSym(
      p = dim(sigma0)[1],
      label = "sigma0"
    )
    for (j in seq_len(dim(sigma0)[2])) {
      for (i in seq_len(dim(sigma0)[1])) {
        if (sigma0[i, j] == 0) {
          params_inicov[i, j] <- "fixed"
        }
      }
    }
  }
  # covariates not allowed at the moment
  return(
    dynr::prep.initial(
      values.inistate = mu0,
      params.inistate = params_inistate,
      values.inicov = sigma0,
      params.inicov = params_inicov
    )
  )
}
