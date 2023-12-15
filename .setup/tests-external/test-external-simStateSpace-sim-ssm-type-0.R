## ---- test-external-simStateSpace-sim-ssm-type-0
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    k <- p <- 3
    iden <- diag(k)
    iden_sqrt <- chol(iden)
    null_vec <- rep(x = 0, times = k)
    mu0 <- null_vec
    sigma0 <- iden
    sigma0_sqrt <- iden_sqrt
    alpha <- null_vec
    beta <- diag(x = 0.50, nrow = k)
    psi <- iden
    psi_sqrt <- iden_sqrt
    nu <- null_vec
    lambda <- iden
    theta <- diag(x = 0.50, nrow = k)
    theta_sqrt <- chol(theta)
    time <- 1000
    burn_in <- 0
    gamma_y <- gamma_eta <- 0.10 * diag(k)
    x <- matrix(
      data = rnorm(n = k * (time + burn_in)),
      ncol = k
    )

    sim <- simStateSpace::SimSSM(
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      type = 0,
      time = time,
      burn_in = burn_in
    )
    data <- simStateSpace::Sim2Matrix(sim)

    obs <- paste0("y", seq_len(k))
    lat <- paste0("eta_", seq_len(p))

    dynr_data <- dynr::dynr.data(
      dataframe = data,
      id = "id",
      time = "time",
      observed = obs
    )

    dynr_initial <- dynr::prep.initial(
      values.inistate = mu0,
      params.inistate = paste0("mu0_", seq_len(k)),
      values.inicov = sigma0,
      params.inicov = matrix(
        data = c(
          "sigma0_11", "fixed",     "fixed",
          "fixed",     "sigma0_22", "fixed",
          "fixed",     "fixed",     "sigma0_33"
        ),
        nrow = 3
      )
    )

    dynr_measurement <- dynr::prep.measurement(
      values.load = diag(k),
      params.load = matrix(data = "fixed", nrow = k, ncol = k),
      state.names = lat,
      obs.names = obs
    )

    dynr_dynamics <- dynr::prep.formulaDynamics(
      formula = list(
        eta_1 ~ beta_11 * eta_1,
        eta_2 ~ beta_22 * eta_2,
        eta_3 ~ beta_33 * eta_3
      ),
      startval = c(
        beta_11 = beta[1, 1],
        beta_22 = beta[2, 2],
        beta_33 = beta[3, 3]
      ),
      isContinuousTime = FALSE
    )

    dynr_noise <- dynr::prep.noise(
      values.latent = psi,
      params.latent = matrix(
        data = c(
          "psi_11", "fixed", "fixed",
          "fixed", "psi_22", "fixed",
          "fixed", "fixed", "psi_33"
        ),
        nrow = 3
      ),
      values.observed = theta,
      params.observed = matrix(
        data = c(
          "theta_11", "fixed", "fixed",
          "fixed", "theta_22", "fixed",
          "fixed", "fixed", "theta_33"
        ),
        nrow = 3
      )
    )

    model <- dynr::dynr.model(
      data = dynr_data,
      initial = dynr_initial,
      measurement = dynr_measurement,
      dynamics = dynr_dynamics,
      noise = dynr_noise
    )
    model@options$maxeval <- 1000000
    model@options$xtol_rel <- .Machine$double.xmin
    model@options$ftol_rel <- .Machine$double.xmin
    model@options$ftol_abs <- .Machine$double.xmin

    model$lb[
      c(
        "beta_11",
        "beta_22",
        "beta_33"
      )
    ] <- -1

    model$ub[
      c(
        "beta_11",
        "beta_22",
        "beta_33"
      )
    ] <- 1

    model$lb[
      c(
        "psi_11",
        "psi_22",
        "psi_33"
      )
    ] <- .Machine$double.xmin

    model$lb[
      c(
        "theta_11",
        "theta_22",
        "theta_33"
      )
    ] <- .Machine$double.xmin

    model$lb[
      c(
        "sigma0_11",
        "sigma0_22",
        "sigma0_33"
      )
    ] <- .Machine$double.xmin

    results <- dynr::dynr.cook(
      model,
      debug_flag = TRUE,
      hessian_flag = FALSE,
      verbose = FALSE
    )

    estimates <- coef(results)
    beta_hat <- matrix(
      data = c(
        estimates[
          c(
            "beta_11",
            "beta_21",
            "beta_31",
            "beta_12",
            "beta_22",
            "beta_32",
            "beta_13",
            "beta_23",
            "beta_33"
          )
        ]
      ),
      nrow = k
    )
  },
  text = "test-external-simStateSpace-sim-ssm-type-0"
)
