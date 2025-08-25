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
        # prepare parameters
        # In this example, the mean vector of the intercept and slope vary.
        # Specifically,
        # there are two sets of values representing two latent classes.
        set.seed(42)
        ## number of individuals
        n <- 10
        ## time points
        time <- 5
        ## dynamic structure
        p <- 2
        mu0_1 <- c(0.615, 1.006) # lower starting point, higher growth
        mu0_2 <- c(1.000, 0.500) # higher starting point, lower growth
        mu0 <- list(mu0_1, mu0_2)
        sigma0 <- matrix(
          data = c(
            1.932,
            0.618,
            0.618,
            0.587
          ),
          nrow = p
        )
        sigma0_l <- list(t(chol(sigma0)))
        ## measurement model
        k <- 1
        theta <- 0.50
        theta_l <- list(sqrt(theta))
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
        ssm <- SimSSMLinGrowthIVary(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
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
        ssm <- SimSSMLinGrowthIVary(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
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
        ssm <- SimSSMLinGrowthIVary(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
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

        testthat::expect_true(
          TRUE
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-lin-growth-i-vary"
)
