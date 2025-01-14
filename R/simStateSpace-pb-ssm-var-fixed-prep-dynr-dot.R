.PBSSMVARFixedPrepDynr <- function(mu0,
                                   sigma0_l,
                                   alpha,
                                   beta,
                                   psi_l,
                                   mu0_fixed,
                                   sigma0_fixed) {
  p <- dim(psi_l)[1]
  k <- p
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
  return(
    list(
      est = est,
      dynr_initial = dynr_initial,
      dynr_measurement = dynr_measurement,
      dynr_noise = dynr_noise,
      dynr_dynamics = dynr_dynamics
    )
  )
}
