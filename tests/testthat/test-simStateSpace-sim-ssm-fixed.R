## ---- test-simStateSpace-sim-ssm-fixed
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
      n = 5,
      time = 100,
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

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)

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
    mu_eta_hat <- colMeans(eta)
    sigma_eta_hat <- stats::cov(eta)
    mu <- simStateSpace:::.Mu0(
      alpha = alpha,
      beta = beta,
      nu = nu
    )
    mu_eta <- as.vector(
      mu$mu_eta
    )
    mu_y <- as.vector(
      mu$mu_y
    )
    sigma <- simStateSpace:::.Sigma0(
      beta = beta,
      psi_l = psi_l,
      lambda = lambda,
      theta_l = theta_l
    )
    sigma_eta <- sigma$sigma_eta
    sigma_y <- sigma$sigma_y
    testthat::test_that(
      paste(text, "mu0", "type 0"),
      {
        testthat::expect_true(
          all(
            abs(
              mu0 - mu0_hat
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "sigma0", "type 0"),
      {
        testthat::expect_true(
          all(
            abs(
              as.vector(
                sigma0 - sigma0_hat
              )
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "mu", "type 0"),
      {
        testthat::expect_true(
          all(
            abs(
              mu_eta - mu_eta_hat
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "sigma", "type 0"),
      {
        testthat::expect_true(
          all(
            abs(
              as.vector(
                sigma_eta - sigma_eta_hat
              )
            ) <= tol
          )
        )
      }
    )
    # Type 1
    ssm <- SimSSMFixed(
      n = 5,
      time = 100,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      alpha = alpha,
      beta = beta,
      psi_l = psi_l,
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
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 2
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
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)
  },
  text = "test-simStateSpace-sim-ssm-fixed",
  tol = 0.01
)
