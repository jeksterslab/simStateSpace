## ---- test-simStateSpace-lin-sde-cov-eta
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
        delta_t <- 0.10
        k <- p <- 3
        iden <- diag(k)
        null_vec <- rep(x = 0, times = k)
        iota <- rep(x = 1, times = p)
        phi <- matrix(
          data = c(
            -0.357,
            0.771,
            -0.450,
            0.0,
            -0.511,
            0.729,
            0,
            0,
            -0.693
          ),
          nrow = k
        )
        sigma <- matrix(
          data = c(
            0.24455556,
            0.02201587,
            -0.05004762,
            0.02201587,
            0.07067800,
            0.01539456,
            -0.05004762,
            0.01539456,
            0.07553061
          ),
          nrow = p
        )
        sigma_l <- t(chol(sigma))
        nu <- rep(x = 1, times = k)
        lambda <- iden
        theta <- iden
        theta_l <- t(chol(theta))
        mu0 <- simStateSpace::LinSDEMeanEta(
          phi = phi,
          iota = iota
        )
        mu <- mu0
        sigma0 <- simStateSpace::LinSDECovEta(
          phi = phi,
          sigma = sigma
        )
        sigma0_l <- t(chol(sigma0))
        sim <- simStateSpace::SimSSMOUFixed(
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
        data <- as.matrix(sim, eta = TRUE)
        eta <- data[, paste0("eta", seq_len(p))]
        testthat::expect_true(
          all(
            (
              c(
                sigma0
              ) - c(
                cov(
                  eta
                )
              )
            ) <= tol
          )
        )
        syl <- simStateSpace:::.SolveSyl(
          A = phi,
          B = t(phi),
          C = sigma
        )
        syl <- (syl + t(syl)) / 2
        testthat::expect_true(
          all(
            (
              c(
                syl
              ) - c(
                sigma0
              )
            ) <= tol
          )
        )
        lya <- simStateSpace:::.SolveLya(
          A = phi,
          Q = sigma
        )
        testthat::expect_true(
          all(
            (
              c(
                lya
              ) - c(
                sigma0
              )
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-cov-eta",
  tol = 0.01
)
