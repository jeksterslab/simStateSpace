## ---- test-simStateSpace-sim-ssm-var
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    k <- 3
    iden <- diag(k)
    iden_sqrt <- chol(iden)
    null_vec <- rep(x = 0, times = k)
    mu0 <- null_vec
    sigma0_sqrt <- iden_sqrt
    alpha <- null_vec
    beta <- diag(x = 0.5, nrow = k)
    psi_sqrt <- iden_sqrt
    time <- 50
    burn_in <- 10
    gamma_eta <- 0.10 * diag(k)
    x <- matrix(
      data = rnorm(n = k * (time + burn_in)),
      ncol = k
    )

    # No covariates
    ssm <- SimSSMVAR(
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
    ssm <- SimSSMVAR(
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
  text = "test-simStateSpace-sim-ssm-var"
)
