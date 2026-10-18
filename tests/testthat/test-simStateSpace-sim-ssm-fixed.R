## ---- test-simStateSpace-sim-ssm-fixed
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        n <- 3
        time <- 8
        delta_t <- 0.5
        p <- 2
        k <- 3
        j <- 2
        mu0 <- c(1, -1)
        sigma0_l <- t(chol(matrix(data = c(1, 0.3, 0.3, 0.5), nrow = p)))
        alpha <- c(0.2, -0.1)
        beta <- matrix(data = c(0.6, 0.1, -0.2, 0.4), nrow = p)
        psi_l <- t(chol(matrix(data = c(0.3, 0.1, 0.1, 0.2), nrow = p)))
        nu <- c(0.5, 1, -0.5)
        lambda <- matrix(data = c(1, 0.8, 0, 0, 0.5, 1), nrow = k)
        theta_l <- diag(c(0.2, 0.3, 0.4))
        set.seed(42)
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            matrix(data = stats::rnorm(n = j * time), nrow = j, ncol = time)
          }
        )
        gamma <- matrix(data = c(0.5, 0, -0.3, 0.2), nrow = p)
        kappa <- matrix(data = c(0.1, 0.2, 0.3, -0.1, 0, 0.4), nrow = k)
        for (type in 0:2) {
          set.seed(type)
          ssm <- SimSSMFixed(
            n = n,
            time = time,
            delta_t = delta_t,
            mu0 = mu0,
            sigma0_l = sigma0_l,
            alpha = alpha,
            beta = beta,
            psi_l = psi_l,
            nu = nu,
            lambda = lambda,
            theta_l = theta_l,
            type = type,
            x = x,
            gamma = gamma,
            kappa = kappa
          )
          set.seed(type)
          noise <- .RefNoise(n = n, time = time, p = p, k = k)
          ref <- .RefSim(
            noise = noise,
            delta_t = delta_t,
            mu0 = mu0,
            sigma0_l = sigma0_l,
            alpha = alpha,
            beta = beta,
            psi_l = psi_l,
            nu = nu,
            lambda = lambda,
            theta_l = theta_l,
            x = if (type > 0) x,
            gamma = if (type > 0) gamma,
            kappa = if (type > 1) kappa
          )
          .RefExpectData(object = ssm, ref = ref)
          .RefExpectMethods(object = ssm, ref = ref)
          .RefExpectModel(
            object = ssm,
            model = "ssm",
            type = type,
            fixed = TRUE,
            fun = "SimSSMFixed"
          )
        }
        testthat::expect_error(
          SimSSMFixed(
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
            type = 3
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-fixed"
)
