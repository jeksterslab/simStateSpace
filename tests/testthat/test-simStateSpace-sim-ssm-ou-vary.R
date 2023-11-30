## ---- test-simStateSpace-sim-ssm-ou-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, phi varies across individuals
    set.seed(42)
    p <- k <- 2
    iden <- diag(p)
    iden_sqrt <- chol(iden)
    n <- 5
    mu0 <- list(c(-3.0, 1.5))
    sigma0_sqrt <- list(iden_sqrt)
    mu <- list(c(5.76, 5.18))
    phi <- list(
      as.matrix(Matrix::expm(diag(x = -0.1, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.2, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.3, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.4, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.5, nrow = k)))
    )
    sigma_sqrt <- list(
      chol(
        matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
      )
    )
    nu <- list(rep(x = 0, times = k))
    lambda <- list(diag(k))
    theta_sqrt <- list(chol(diag(x = 0.50, nrow = k)))
    delta_t <- 0.10
    time <- 50
    burn_in <- 10
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
    ssm <- SimSSMOUVary(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      mu = mu,
      phi = phi,
      sigma_sqrt = sigma_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      type = 0,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Type 1
    ssm <- SimSSMOUVary(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      mu = mu,
      phi = phi,
      sigma_sqrt = sigma_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      gamma_eta = gamma_eta,
      x = x,
      type = 1,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Type 2
    ssm <- SimSSMOUVary(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      mu = mu,
      phi = phi,
      sigma_sqrt = sigma_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = 2,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

  },
  text = "test-simStateSpace-sim-ssm-ou-vary"
)
