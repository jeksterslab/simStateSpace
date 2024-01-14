## ---- test-simStateSpace-sim-ssm-lin-growth-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, the mean vector of the intercept and slope vary.
    # Specifically,
    # there are two sets of values representing two latent classes.
    set.seed(42)
    n <- 10
    mu0_1 <- c(0.615, 1.006) # lower starting point, higher growth
    mu0_2 <- c(1.000, 0.500) # higher starting point, lower growth
    mu0 <- list(mu0_1, mu0_2)
    sigma0 <- list(
      matrix(
        data = c(
          1.932,
          0.618,
          0.618,
          0.587
        ),
        nrow = 2
      )
    )
    theta <- list(0.6)
    time <- 10
    gamma_y <- list(matrix(data = 0.10, nrow = 1, ncol = 2))
    gamma_eta <- list(matrix(data = 0.10, nrow = 2, ncol = 2))
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        return(
          matrix(
            data = rnorm(n = 2 * time),
            ncol = 2
          )
        )
      }
    )

    # Type 0
    ssm <- SimSSMLinGrowthIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      theta = theta,
      type = 0,
      time = time
    )

    as.data.frame(ssm, eta = TRUE)
    as.data.frame(ssm, eta = FALSE)
    as.data.frame(ssm, eta = TRUE, long = FALSE)
    as.data.frame(ssm, eta = FALSE, long = FALSE)
    as.matrix(ssm, eta = TRUE)
    as.matrix(ssm, eta = FALSE)
    as.matrix(ssm, eta = TRUE, long = FALSE)
    as.matrix(ssm, eta = FALSE, long = FALSE)
    print(ssm)
    plot(ssm)
    plot(ssm, eta = TRUE)

    # Type 1
    ssm <- SimSSMLinGrowthIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      theta = theta,
      gamma_eta = gamma_eta,
      x = x,
      type = 1,
      time = time
    )

    as.data.frame(ssm, eta = TRUE)
    as.data.frame(ssm, eta = FALSE)
    as.data.frame(ssm, eta = TRUE, long = FALSE)
    as.data.frame(ssm, eta = FALSE, long = FALSE)
    as.matrix(ssm, eta = TRUE)
    as.matrix(ssm, eta = FALSE)
    as.matrix(ssm, eta = TRUE, long = FALSE)
    as.matrix(ssm, eta = FALSE, long = FALSE)
    print(ssm)
    plot(ssm)
    plot(ssm, eta = TRUE)

    # Type 2
    ssm <- SimSSMLinGrowthIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      theta = theta,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = 2,
      time = time
    )

    as.data.frame(ssm, eta = TRUE)
    as.data.frame(ssm, eta = FALSE)
    as.data.frame(ssm, eta = TRUE, long = FALSE)
    as.data.frame(ssm, eta = FALSE, long = FALSE)
    as.matrix(ssm, eta = TRUE)
    as.matrix(ssm, eta = FALSE)
    as.matrix(ssm, eta = TRUE, long = FALSE)
    as.matrix(ssm, eta = FALSE, long = FALSE)
    print(ssm)
    plot(ssm)
    plot(ssm, eta = TRUE)

    # Error
    testthat::test_that(
      paste(text, "error"),
      {
        testthat::expect_error(
          SimSSMLinGrowthIVary(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            theta = theta,
            gamma_y = gamma_y,
            gamma_eta = gamma_eta,
            x = x,
            type = 3,
            time = time
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-lin-growth-i-vary"
)
