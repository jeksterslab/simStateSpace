## ---- test-simStateSpace-sim-ssm-ou-fixed
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    p <- k <- 2
    iden <- diag(p)
    iden_sqrt <- chol(iden)
    n <- 5
    mu0 <- c(-3.0, 1.5)
    sigma0_sqrt <- iden_sqrt
    mu <- c(5.76, 5.18)
    phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
    sigma_sqrt <- chol(
      matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
    )
    nu <- rep(x = 0, times = k)
    lambda <- diag(k)
    theta_sqrt <- chol(diag(x = 0.50, nrow = k))
    delta_t <- 0.10
    time <- 50
    burn_in <- 10

    ssm <- SimSSMOUFixed(
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

    ssm <- SimSSMOU(
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
  text = "test-simStateSpace-sim-ssm-ou-fixed"
)
