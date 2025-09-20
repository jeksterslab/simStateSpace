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
        # prepare parameters
        # In this example, beta varies across individuals.
        set.seed(42)
        ## number of individuals
        n <- 5
        ## time points
        time <- 50
        ## dynamic structure
        p <- 3
        mu0 <- list(
          rep(x = 0, times = p)
        )
        sigma0 <- 0.001 * diag(p)
        sigma0_l <- list(
          t(chol(sigma0))
        )
        alpha <- list(
          rep(x = 0, times = p)
        )
        beta <- SimBetaN(
          n = 5,
          beta = 0.1 * diag(p),
          vcov_beta_vec_l = t(
            chol(
              0.1 * diag(p^2)
            )
          )
        )
        psi <- 0.001 * diag(p)
        psi_l <- list(
          t(chol(psi))
        )
        ## measurement model
        k <- 3
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
        ssm <- SimSSMIVary(
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
        ssm <- SimSSMIVary(
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
        ssm <- SimSSMIVary(
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
  text = "test-simStateSpace-sim-ssm-i-vary"
)
