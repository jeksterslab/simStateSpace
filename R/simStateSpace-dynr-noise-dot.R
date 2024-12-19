.DynrNoise <- function(process_l,
                       theta_l,
                       continuous) {
  .GenerateLabels <- function(x,
                              label) {
    params <- .LabelSym(
      p = dim(x)[1],
      label = label
    )
    zero_indices <- which(
      x = x == 0,
      arr.ind = TRUE
    )
    params[zero_indices] <- "fixed"
    return(params)
  }
  theta <- tcrossprod(theta_l)
  process <- tcrossprod(process_l)
  params_observed <- .GenerateLabels(
    x = theta,
    label = "theta"
  )
  params_latent <- .GenerateLabels(
    x = process,
    label = if (continuous) "sigma" else "psi"
  )
  return(
    dynr::prep.noise(
      values.latent = process,
      params.latent = params_latent,
      values.observed = theta,
      params.observed = params_observed
    )
  )
}
