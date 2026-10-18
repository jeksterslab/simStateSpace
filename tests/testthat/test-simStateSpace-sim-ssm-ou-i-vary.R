## ---- test-simStateSpace-sim-ssm-ou-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # parameters are lists of length 2, recycled over n = 5;
        # the second person has mu = 0 and sigma_l = 0
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
        mu <- list(c(0.5, -0.3), c(0, 0))
        phi <- list(
          # complex eigenvalues -0.4 +/- 0.6i
          matrix(data = c(-0.4, -0.6, 0.6, -0.4), nrow = p),
          # real eigenvalues
          matrix(data = c(-0.5, 0.2, 0.1, -0.8), nrow = p)
        )
        sigma_l <- list(
          t(chol(matrix(data = c(0.3, 0.1, 0.1, 0.2), nrow = p))),
          matrix(data = 0, nrow = p, ncol = p)
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
        disc <- lapply(
          X = 1:2,
          FUN = function(i) {
            .RefDiscretize(
              iota = -phi[[i]] %*% mu[[i]],
              phi = phi[[i]],
              sigma = sigma_l[[i]] %*% t(sigma_l[[i]]),
              delta_t = delta_t
            )
          }
        )
        for (type in 0:2) {
          set.seed(type)
          ssm <- SimSSMOUIVary(
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
            alpha = lapply(X = disc, FUN = "[[", "alpha"),
            beta = lapply(X = disc, FUN = "[[", "beta"),
            psi_l = lapply(X = disc, FUN = "[[", "psi_l"),
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
            model = "ou",
            type = type,
            fixed = FALSE,
            fun = "SimSSMOUIVary"
          )
        }
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-ou-i-vary"
)
