## ---- test-simStateSpace-sim-ssm-lin-growth-i-vary
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
        j <- 2
        mu0 <- list(c(0.6, 0.3), c(-0.2, 0.1))
        sigma0_l <- list(
          t(chol(matrix(data = c(1, 0.2, 0.2, 0.5), nrow = 2))),
          diag(c(0.3, 0.2))
        )
        theta_l <- list(0.5, 0.2)
        set.seed(42)
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            matrix(data = stats::rnorm(n = j * time), nrow = j, ncol = time)
          }
        )
        gamma <- list(
          matrix(data = c(0.5, 0, -0.3, 0.2), nrow = 2),
          matrix(data = c(0, 0.4, 0.1, 0), nrow = 2)
        )
        kappa <- list(
          matrix(data = c(0.1, -0.2), nrow = 1),
          matrix(data = c(0.3, 0), nrow = 1)
        )
        for (type in 0:2) {
          set.seed(type)
          ssm <- SimSSMLinGrowthIVary(
            n = n,
            time = time,
            mu0 = mu0,
            sigma0_l = sigma0_l,
            theta_l = theta_l,
            type = type,
            x = x,
            gamma = gamma,
            kappa = kappa
          )
          set.seed(type)
          noise <- .RefNoise(n = n, time = time, p = 2, k = 1)
          ref <- .RefSim(
            noise = noise,
            delta_t = 1,
            mu0 = mu0,
            sigma0_l = sigma0_l,
            alpha = c(0, 0),
            beta = rbind(c(1, 1), c(0, 1)),
            psi_l = matrix(data = 0, nrow = 2, ncol = 2),
            nu = 0,
            lambda = matrix(data = c(1, 0), nrow = 1),
            theta_l = lapply(X = theta_l, FUN = as.matrix),
            x = if (type > 0) x,
            gamma = if (type > 0) gamma,
            kappa = if (type > 1) kappa
          )
          .RefExpectData(object = ssm, ref = ref)
          .RefExpectMethods(object = ssm, ref = ref)
          .RefExpectModel(
            object = ssm,
            model = "lingrowth",
            type = type,
            fixed = FALSE,
            fun = "SimSSMLinGrowthIVary"
          )
        }
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-lin-growth-i-vary"
)
