## ---- test-simStateSpace-sim-ssm-lin-sde-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # prepare parameters
        # In this example, phi varies across individuals.
        set.seed(42)
        ## number of individuals
        n <- 5
        ## time points
        time <- 50
        delta_t <- 0.10
        ## dynamic structure
        p <- 2
        mu0 <- list(
          c(-3.0, 1.5)
        )
        sigma0 <- 0.001 * diag(p)
        sigma0_l <- list(
          t(chol(sigma0))
        )
        iota <- list(
          c(0.317, 0.230)
        )
        phi <- SimPhiN(
          n = 5,
          phi = -0.1 * diag(p),
          vcov_phi_vec_l = t(
            chol(
              0.1 * diag(p^2)
            )
          )
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
        sigma_l <- list(
          t(chol(sigma))
        )
        ## measurement model
        k <- 2
        nu <- list(
          rep(x = 0, times = k)
        )
        lambda <- list(
          diag(k)
        )
        theta <- 0.001 * diag(k)
        theta_l <- list(
          t(chol(theta))
        )
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
        gamma <- list(
          diag(x = 0.10, nrow = p, ncol = j)
        )
        kappa <- list(
          diag(x = 0.10, nrow = k, ncol = j)
        )

        # Type 0
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
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
        as.data.frame.simstatespace(ssm, burnin = 5, reset_time = TRUE)
        as.data.frame.simstatespace(ssm, burnin = 5, reset_time = FALSE)
        as.matrix.simstatespace(ssm, eta = TRUE)
        as.matrix.simstatespace(ssm, eta = FALSE)
        as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
        as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
        as.matrix.simstatespace(ssm, burnin = 5, reset_time = TRUE)
        as.matrix.simstatespace(ssm, burnin = 5, reset_time = FALSE)
        print.simstatespace(ssm)
        plot.simstatespace(ssm, id = 1:3, time = 0:4)
        plot.simstatespace(ssm, eta = TRUE)
        plot.simstatespace(ssm, burnin = 5, reset_time = FALSE)

        # Type 1
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
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
        as.data.frame.simstatespace(ssm, burnin = 5, reset_time = TRUE)
        as.data.frame.simstatespace(ssm, burnin = 5, reset_time = FALSE)
        as.matrix.simstatespace(ssm, eta = TRUE)
        as.matrix.simstatespace(ssm, eta = FALSE)
        as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
        as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
        as.matrix.simstatespace(ssm, burnin = 5, reset_time = TRUE)
        as.matrix.simstatespace(ssm, burnin = 5, reset_time = FALSE)
        print.simstatespace(ssm)
        plot.simstatespace(ssm, id = 1:3, time = 0:4)
        plot.simstatespace(ssm, eta = TRUE)
        plot.simstatespace(ssm, burnin = 5, reset_time = FALSE)

        # Type 2
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
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
        as.data.frame.simstatespace(ssm, burnin = 5, reset_time = TRUE)
        as.data.frame.simstatespace(ssm, burnin = 5, reset_time = FALSE)
        as.matrix.simstatespace(ssm, eta = TRUE)
        as.matrix.simstatespace(ssm, eta = FALSE)
        as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
        as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
        as.matrix.simstatespace(ssm, burnin = 5, reset_time = TRUE)
        as.matrix.simstatespace(ssm, burnin = 5, reset_time = FALSE)
        print.simstatespace(ssm)
        plot.simstatespace(ssm, id = 1:3, time = 0:4)
        plot.simstatespace(ssm, eta = TRUE)
        plot.simstatespace(ssm, burnin = 5, reset_time = FALSE)

        # sigma_l zero --------------------------------------------
        sigma_l <- list(
          matrix(
            data = 0,
            nrow = p,
            ncol = p
          )
        )

        # Type 0
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
          phi = phi,
          sigma_l = sigma_l,
          nu = nu,
          lambda = lambda,
          theta_l = theta_l,
          type = 0
        )

        # Type 1
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
          phi = phi,
          sigma_l = sigma_l,
          nu = nu,
          lambda = lambda,
          theta_l = theta_l,
          type = 1,
          x = x,
          gamma = gamma
        )

        # Type 2
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
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

        testthat::expect_error(
          as.matrix.simstatespace(ssm, burnin = time),
          regexp = paste0(
            "`burnin` should not be greater than the measurement occasions.\n"
          )
        )

        testthat::expect_true(
          TRUE
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-lin-sde-i-vary"
)
