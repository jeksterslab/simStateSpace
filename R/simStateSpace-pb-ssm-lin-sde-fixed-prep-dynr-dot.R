.PBSSMLinSDEFixedPrepDynr <- function(mu0,
                                      sigma0_l,
                                      iota,
                                      phi,
                                      sigma_l,
                                      nu,
                                      lambda,
                                      theta_l,
                                      mu0_fixed,
                                      sigma0_fixed) {
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
    lambda = lambda,
    nu = nu
  )
  nu_values <- .Vec(
    dynr_measurement$values.int[[1]]
  )
  nu_labels <- .Vec(
    dynr_measurement$params.int[[1]]
  )
  names(nu_values) <- nu_labels
  dynr_noise <- .DynrNoise(
    process_l = sigma_l,
    theta_l = theta_l,
    continuous = TRUE
  )
  sigma_values <- .Vec(
    dynr_noise$values.latent[[1]]
  )
  sigma_labels <- .Vec(
    dynr_noise$params.latent[[1]]
  )
  names(sigma_values) <- sigma_labels
  theta_values <- .Vech(
    dynr_noise$values.observed[[1]]
  )
  theta_labels <- .Vech(
    dynr_noise$params.observed[[1]]
  )
  names(theta_values) <- theta_labels
  dynr_dynamics <- .DynrDynamics(
    dynamics = phi,
    intercept = iota,
    continuous = TRUE,
    ou = FALSE
  )
  dynamics_values <- dynr_dynamics$startval
  dynr_dynamics <- dynr_dynamics$dynamics
  est <- c(
    dynamics_values,
    sigma_values,
    nu_values,
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
