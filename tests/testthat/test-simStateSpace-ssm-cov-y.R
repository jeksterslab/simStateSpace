## ---- test-simStateSpace-ssm-cov-y
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        set.seed(42)
        n <- 1000
        time <- 1000
        k <- p <- 3
        iden <- diag(k)
        null_vec <- rep(x = 0, times = k)
        alpha <- null_vec
        beta <- matrix(
          data = c(
            0.7,
            0.5,
            -0.1,
            0.0,
            0.6,
            0.4,
            0,
            0,
            0.5
          ),
          nrow = k
        )
        psi <- 0.1 * iden
        psi_l <- t(chol(psi))
        nu <- rep(x = 1, times = k)
        theta <- iden
        theta <- 0.2 * iden
        theta_l <- t(chol(theta))
        mu0 <- simStateSpace::SSMMeanEta(
          beta = beta,
          alpha = alpha
        )
        sigma0 <- simStateSpace::SSMCovEta(
          beta = beta,
          psi = psi
        )
        sigma0_l <- t(chol(sigma0))
        sigma_y <- simStateSpace::SSMCovY(
          lambda = lambda,
          theta = theta,
          cov_eta = sigma0
        )
        sim <- simStateSpace::SimSSMFixed(
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
        data <- as.matrix(sim, eta = TRUE)
        y <- data[, paste0("y", seq_len(p))]
        testthat::expect_true(
          all(
            (
              c(
                sigma0
              ) - c(
                cov(
                  y
                )
              )
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-ssm-cov-y",
  tol = 0.01
)
