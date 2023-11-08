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
    n <- 3
    mu0 <- list(c(-3.0, 1.5))
    sigma0_sqrt <- list(iden_sqrt)
    mu <- list(c(5.76, 5.18))
    phi <- list(
      as.matrix(Matrix::expm(diag(x = -0.1, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.2, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.3, nrow = k)))
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
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
  },
  text = "test-simStateSpace-sim-ssm-ou-vary"
)
