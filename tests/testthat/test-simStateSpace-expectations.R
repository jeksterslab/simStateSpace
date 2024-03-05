## ---- test-simStateSpace-expectations
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    # prepare parameters
    set.seed(42)
    ## number of individuals
    n <- 1000
    ## time points
    time <- 1000
    ## dynamic structure
    p <- 3
    mu0 <- rep(x = 0, times = p)
    sigma0 <- 0.001 * diag(p)
    sigma0_l <- t(chol(sigma0))
    alpha <- rep(x = 0, times = p)
    beta <- 0.50 * diag(p)
    psi <- 0.001 * diag(p)
    psi_l <- t(chol(psi))
    ## measurement model
    k <- 3
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
    ssm <- SimSSMFixed(
      n = n,
      time = time,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      alpha = alpha,
      beta = beta,
      psi_l = psi_l,
      nu = nu,
      lambda = lambda,
      theta_l = theta_l,
      type = 0
    )
    data <- as.matrix.simstatespace(ssm, eta = TRUE)
    eta <- data[, c("time", paste0("eta", 1:p))]
    eta0 <- eta[which(eta[, "time"] == 0), ]
    eta0 <- eta0[, paste0("eta", 1:p)]
    eta <- eta[, paste0("eta", 1:p)]

    mu0_hat <- colMeans(eta0)
    sigma0_hat <- stats::cov(eta0)
    mu_hat <- colMeans(eta)
    sigma_hat <- stats::cov(eta)

    mu <- MuEta0(beta = beta, alpha = alpha)
    sigma <- SigmaEta0(beta = beta, psi_l = psi_l)
  },
  text = "test-simStateSpace-expectations",
  tol = 0.001
)
