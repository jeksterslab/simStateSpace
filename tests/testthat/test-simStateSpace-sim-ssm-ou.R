## ---- test-simStateSpace-sim-ssm-ou
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
    gamma_y <- gamma_eta <- 0.10 * diag(k)
    x <- matrix(
      data = rnorm(n = k * (time + burn_in)),
      ncol = k
    )

    # Type 0
    ssm <- SimSSMOU(
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
    ssm <- SimSSMOU(
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
    ssm <- SimSSMOU(
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
  text = "test-simStateSpace-sim-ssm-ou"
)
