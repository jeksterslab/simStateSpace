## ---- test-simStateSpace-sim-ssm-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # parameters are lists of length 2, recycled over n = 5
        n <- 5
        time <- 8
        delta_t <- 0.5
        p <- 2
        k <- 3
        j <- 2
        mu0 <- list(c(1, -1), c(0, 0.5))
        sigma0_l <- list(
          t(chol(matrix(data = c(1, 0.3, 0.3, 0.5), nrow = p))),
          diag(c(0.5, 0.2))
        )
        alpha <- list(c(0.2, -0.1), c(0, 0.3))
        beta <- list(
          matrix(data = c(0.6, 0.1, -0.2, 0.4), nrow = p),
          matrix(data = c(0.3, -0.2, 0.1, 0.5), nrow = p)
        )
        psi_l <- list(
          t(chol(matrix(data = c(0.3, 0.1, 0.1, 0.2), nrow = p))),
          diag(c(0.1, 0.4))
        )
        nu <- list(c(0.5, 1, -0.5), c(0, 0, 1))
        lambda <- list(
          matrix(data = c(1, 0.8, 0, 0, 0.5, 1), nrow = k),
          matrix(data = c(1, 0, 0.3, 0, 1, 0.7), nrow = k)
        )
        theta_l <- list(diag(c(0.2, 0.3, 0.4)), diag(c(0.5, 0.1, 0.2)))
        set.seed(42)
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            matrix(data = stats::rnorm(n = j * time), nrow = j, ncol = time)
          }
        )
        gamma <- list(
          matrix(data = c(0.5, 0, -0.3, 0.2), nrow = p),
          matrix(data = c(0, 0.4, 0.1, 0), nrow = p)
        )
        kappa <- list(
          matrix(data = c(0.1, 0.2, 0.3, -0.1, 0, 0.4), nrow = k),
          matrix(data = c(0, 0.3, 0, 0.2, 0.1, 0), nrow = k)
        )
        for (type in 0:2) {
          set.seed(type)
          ssm <- SimSSMIVary(
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
            fixed = FALSE,
            fun = "SimSSMIVary"
          )
        }
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-i-vary"
)
