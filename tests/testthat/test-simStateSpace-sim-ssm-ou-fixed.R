## ---- test-simStateSpace-sim-ssm-ou-fixed
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    p <- k <- 2
    iden <- diag(p)
    n <- 5
    mu0 <- c(-3.0, 1.5)
    sigma0 <- iden
    mu <- c(5.76, 5.18)
    phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
    sigma <- matrix(
      data = c(2.79, 0.06, 0.06, 3.27),
      nrow = p
    )
    nu <- rep(x = 0, times = k)
    lambda <- diag(k)
    theta <- diag(x = 0.50, nrow = k)
    delta_t <- 0.10
    time <- 50
    burn_in <- 10
    gamma_y <- gamma_eta <- 0.10 * diag(k)
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        return(
          matrix(
            data = rnorm(n = k * (time + burn_in)),
            ncol = k
          )
        )
      }
    )

    # Type 0
    ssm <- simStateSpace::SimSSMOUFixed(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      mu = mu,
      phi = phi,
      sigma = sigma,
      nu = nu,
      lambda = lambda,
      theta = theta,
      type = 0,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = (0:4) * 0.10)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 1
    ssm <- simStateSpace::SimSSMOUFixed(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      mu = mu,
      phi = phi,
      sigma = sigma,
      nu = nu,
      lambda = lambda,
      theta = theta,
      gamma_eta = gamma_eta,
      x = x,
      type = 1,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = (0:4) * 0.10)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 2
    ssm <- simStateSpace::SimSSMOUFixed(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      mu = mu,
      phi = phi,
      sigma = sigma,
      nu = nu,
      lambda = lambda,
      theta = theta,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = 2,
      delta_t = delta_t,
      time = time,
      burn_in = burn_in
    )

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = (0:4) * 0.10)
    plot.simstatespace(ssm, eta = TRUE)

    # Error
    testthat::test_that(
      paste(text, "error"),
      {
        testthat::expect_error(
          simStateSpace::SimSSMOUFixed(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            mu = mu,
            phi = phi,
            sigma = sigma,
            nu = nu,
            lambda = lambda,
            theta = theta,
            gamma_y = gamma_y,
            gamma_eta = gamma_eta,
            x = x,
            type = 3,
            delta_t = delta_t,
            time = time,
            burn_in = burn_in
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "error type 1"),
      {
        testthat::expect_error(
          simStateSpace::SimSSMOUFixed(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            mu = mu,
            phi = phi,
            sigma = sigma,
            nu = nu,
            lambda = lambda,
            theta = theta,
            type = 1,
            delta_t = delta_t,
            time = time,
            burn_in = burn_in
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "error type 2"),
      {
        testthat::expect_error(
          simStateSpace::SimSSMOUFixed(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            mu = mu,
            phi = phi,
            sigma = sigma,
            nu = nu,
            lambda = lambda,
            theta = theta,
            type = 2,
            delta_t = delta_t,
            time = time,
            burn_in = burn_in
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-ou-fixed"
)
