## ---- test-simStateSpace-sim-ssm-ou-fixed
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    set.seed(42)
    ## number of individuals
    n <- 5
    ## time points
    time <- 50
    delta_t <- 0.10
    ## dynamic structure
    p <- 2
    mu0 <- c(-3.0, 1.5)
    sigma0 <- 0.001 * diag(p)
    sigma0_l <- t(chol(sigma0))
    mu <- c(5.76, 5.18)
    phi <- matrix(
      data = c(
        -0.10,
        0.05,
        0.05,
        -0.10
      ),
      nrow = p
    )
    sigma <- matrix(
      data = c(
        2.79,
        0.06,
        0.06,
        3.27
      ),
      nrow = p
    )
    sigma_l <- t(chol(sigma))
    ## measurement model
    k <- 2
    nu <- rep(x = 0, times = k)
    lambda <- diag(k)
    theta <- 0.001 * diag(k)
    theta_l <- t(chol(theta))
    ## covariates
    j <- 2
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        matrix(
          data = stats::rnorm(n = time * j),
          nrow = j,
          ncol = time
        )
      }
    )
    gamma <- diag(x = 0.10, nrow = p, ncol = j)
    kappa <- diag(x = 0.10, nrow = k, ncol = j)

    # Type 0
    ssm <- SimSSMOUFixed(
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      mu = mu,
      phi = phi,
      sigma_l = sigma_l,
      nu = nu,
      lambda = lambda,
      theta_l = theta_l,
      type = 0
    )

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = (0:4) * 0.10)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 1
    ssm <- SimSSMOUFixed(
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      mu = mu,
      phi = phi,
      sigma_l = sigma_l,
      nu = nu,
      lambda = lambda,
      theta_l = theta_l,
      type = 1,
      x = x,
      gamma = gamma
    )

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = (0:4) * 0.10)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 2
    ssm <- SimSSMOUFixed(
      n = n,
      time = time,
      delta_t = delta_t,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      mu = mu,
      phi = phi,
      sigma_l = sigma_l,
      nu = nu,
      lambda = lambda,
      theta_l = theta_l,
      type = 2,
      x = x,
      gamma = gamma,
      kappa = kappa
    )

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = (0:4) * 0.10)
    plot.simstatespace(ssm, eta = TRUE)
  },
  text = "test-simStateSpace-sim-ssm-ou-fixed"
)
