## ---- test-simStateSpace-sim-ssm-var-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, beta varies across individuals
    set.seed(42)
    k <- 3
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
    time <- 50
    burn_in <- 10
    gamma_eta <- list(0.10 * diag(k))
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

    # No covariates
    ssm <- SimSSMVARVary(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # With covariates
    ssm <- SimSSMVARVary(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      gamma_eta = gamma_eta,
      x = x,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

  },
  text = "test-simStateSpace-sim-ssm-var-vary"
)
