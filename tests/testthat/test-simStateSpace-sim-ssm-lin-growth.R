## ---- test-simStateSpace-sim-ssm-lin-growth
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    n <- 5
    mu0 <- c(0.615, 1.006)
    sigma0 <- matrix(
      data = c(
        1.932,
        0.618,
        0.618,
        0.587
      ),
      nrow = 2
    )
    sigma0_sqrt <- chol(sigma0)
    theta <- 0.6
    theta_sqrt <- sqrt(theta)
    time <- 10
    gamma_y <- matrix(data = 0.10, nrow = 1, ncol = 2)
    gamma_eta <- matrix(data = 0.10, nrow = 2, ncol = 2)
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        return(
          matrix(
            data = rnorm(n = 2 * time),
            ncol = 2
          )
        )
      }
    )

    # Type 0
    ssm <- SimSSMLinGrowth(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      theta_sqrt = theta_sqrt,
      type = 0,
      time = time
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Type 1
    ssm <- SimSSMLinGrowth(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      theta_sqrt = theta_sqrt,
      gamma_eta = gamma_eta,
      x = x,
      type = 1,
      time = time
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Type 2
    ssm <- SimSSMLinGrowth(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      theta_sqrt = theta_sqrt,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = 2,
      time = time
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)
  },
  text = "test-simStateSpace-sim-ssm-lin-growth"
)
