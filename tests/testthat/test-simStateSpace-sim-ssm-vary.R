## ---- test-simStateSpace-sim-ssm-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, beta varies across individuals
    set.seed(42)
    k <- p <- 3
    iden <- diag(k)
    iden_sqrt <- chol(iden)
    null_vec <- rep(x = 0, times = k)
    n <- 5
    mu0 <- list(null_vec)
    sigma0_sqrt <- list(iden_sqrt)
    alpha <- list(null_vec)
    beta <- list(
      diag(x = 0.1, nrow = k),
      diag(x = 0.2, nrow = k),
      diag(x = 0.3, nrow = k),
      diag(x = 0.4, nrow = k),
      diag(x = 0.5, nrow = k)
    )
    psi_sqrt <- list(iden_sqrt)
    nu <- list(null_vec)
    lambda <- list(iden)
    theta_sqrt <- list(chol(diag(x = 0.50, nrow = k)))
    time <- 50
    burn_in <- 0
    gamma_y <- gamma_eta <- list(0.10 * diag(k))
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        return(
          matrix(
            data = rnorm(n = k * (time + burn_in)),
            ncol = k
          )
        )
      }
    )

    # Type 0
    ssm <- SimSSMVary(
      n = n,
      type = 0,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)

    # Type 1
    ssm <- SimSSMVary(
      n = n,
      type = 1,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      gamma_eta = gamma_eta,
      x = x,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)

    # Type 2
    ssm <- SimSSMVary(
      n = n,
      type = 2,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
  },
  text = "test-simStateSpace-sim-ssm-vary"
)
