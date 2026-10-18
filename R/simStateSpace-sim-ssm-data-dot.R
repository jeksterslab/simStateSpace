.SimSSMData <- function(n, time, delta_t,
                        mu0, sigma0_l,
                        alpha, beta, psi_l,
                        nu = NULL, lambda = NULL, theta_l = NULL,
                        type = 0,
                        x = NULL, gamma = NULL, kappa = NULL,
                        fixed = TRUE) {
  # Simulate the data of every state space model with one kernel.
  # Parameters are lists of parameter sets; person i uses set
  # i modulo the length of the list. With `fixed = TRUE`, each parameter
  # is one set. `x` is a list of covariate matrices in both cases.
  # Leave `nu`, `lambda`, and `theta_l` as `NULL` for latent-only models.
  # The callers check the inputs first (`.CheckSimulation()`).
  if (type < 1) {
    gamma <- NULL
  }
  if (type < 2) {
    kappa <- NULL
  }
  par <- list(
    mu0 = mu0, sigma0_l = sigma0_l,
    alpha = alpha, beta = beta, psi_l = psi_l,
    nu = nu, lambda = lambda, theta_l = theta_l,
    gamma = gamma, kappa = kappa
  )
  par <- lapply(
    X = par,
    FUN = function(i) {
      if (is.null(i)) {
        list()
      } else if (fixed) {
        list(i)
      } else {
        i
      }
    }
  )
  .SimSSMKernel(
    n = n,
    time = time,
    delta_t = delta_t,
    mu0 = par$mu0,
    sigma0_l = par$sigma0_l,
    alpha = par$alpha,
    beta = par$beta,
    psi_l = par$psi_l,
    nu = par$nu,
    lambda = par$lambda,
    theta_l = par$theta_l,
    x = if (type > 0) x else list(),
    gamma = par$gamma,
    kappa = par$kappa
  )
}
