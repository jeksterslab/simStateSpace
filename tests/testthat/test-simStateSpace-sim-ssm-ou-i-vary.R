## ---- test-simStateSpace-sim-ssm-ou-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, phi varies across individuals
    set.seed(42)
    p <- k <- 2
    iden <- diag(p)
    n <- 5
    mu0 <- list(c(-3.0, 1.5))
    sigma0 <- list(iden)
    mu <- list(c(5.76, 5.18))
    phi <- list(
      as.matrix(Matrix::expm(diag(x = -0.1, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.2, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.3, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.4, nrow = k))),
      as.matrix(Matrix::expm(diag(x = -0.5, nrow = k)))
    )
    sigma <- list(
      matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
    )
    nu <- list(rep(x = 0, times = k))
    lambda <- list(diag(k))
    theta <- list(diag(x = 0.50, nrow = k))
    delta_t <- 0.10
    time <- 50
    burn_in <- 10
    gamma_y <- gamma_eta <- list(0.10 * diag(k))
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
    ssm <- SimSSMOUIVary(
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

    as.data.frame(ssm, eta = TRUE)
    as.data.frame(ssm, eta = FALSE)
    as.data.frame(ssm, eta = TRUE, long = FALSE)
    as.data.frame(ssm, eta = FALSE, long = FALSE)
    as.matrix(ssm, eta = TRUE)
    as.matrix(ssm, eta = FALSE)
    as.matrix(ssm, eta = TRUE, long = FALSE)
    as.matrix(ssm, eta = FALSE, long = FALSE)
    print(ssm)

    # Type 1
    ssm <- SimSSMOUIVary(
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

    as.data.frame(ssm, eta = TRUE)
    as.data.frame(ssm, eta = FALSE)
    as.data.frame(ssm, eta = TRUE, long = FALSE)
    as.data.frame(ssm, eta = FALSE, long = FALSE)
    as.matrix(ssm, eta = TRUE)
    as.matrix(ssm, eta = FALSE)
    as.matrix(ssm, eta = TRUE, long = FALSE)
    as.matrix(ssm, eta = FALSE, long = FALSE)
    print(ssm)

    # Type 2
    ssm <- SimSSMOUIVary(
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

    as.data.frame(ssm, eta = TRUE)
    as.data.frame(ssm, eta = FALSE)
    as.data.frame(ssm, eta = TRUE, long = FALSE)
    as.data.frame(ssm, eta = FALSE, long = FALSE)
    as.matrix(ssm, eta = TRUE)
    as.matrix(ssm, eta = FALSE)
    as.matrix(ssm, eta = TRUE, long = FALSE)
    as.matrix(ssm, eta = FALSE, long = FALSE)
    print(ssm)

    # Error
    testthat::test_that(
      paste(text, "error"),
      {
        testthat::expect_error(
          SimSSMOUIVary(
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
  },
  text = "test-simStateSpace-sim-ssm-ou-i-vary"
)