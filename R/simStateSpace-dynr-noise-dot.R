.DynrNoise <- function(process_l,
                       theta_l,
                       continuous) {
  theta <- tcrossprod(theta_l)
  params_observed <- .LabelSym(
    p = dim(theta)[1],
    label = "theta"
  )
  for (j in seq_len(dim(theta)[2])) {
    for (i in seq_len(dim(theta)[1])) {
      if (theta[i, j] == 0) {
        params_observed[i, j] <- "fixed"
      }
    }
  }
  process <- tcrossprod(process_l)
  if (continuous) {
    params_latent <- .LabelSym(
      p = dim(process)[1],
      label = "sigma"
    )
  } else {
    params_latent <- .LabelSym(
      p = dim(process)[1],
      label = "psi"
    )
  }
  for (j in seq_len(dim(process)[2])) {
    for (i in seq_len(dim(process)[1])) {
      if (process[i, j] == 0) {
        params_latent[i, j] <- "fixed"
      }
    }
  }
  return(
    dynr::prep.noise(
      values.latent = process,
      params.latent = params_latent,
      values.observed = theta,
      params.observed = params_observed
    )
  )
}
