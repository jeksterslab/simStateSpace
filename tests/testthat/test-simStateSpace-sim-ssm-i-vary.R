## ---- test-simStateSpace-sim-ssm-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, beta varies across individuals
    set.seed(42)
    k <- p <- 3
    iden <- diag(k)
    null_vec <- rep(x = 0, times = k)
    n <- 5
    mu0 <- list(null_vec)
    sigma0 <- list(iden)
    alpha <- list(null_vec)
    beta <- list(
      diag(x = 0.1, nrow = k),
      diag(x = 0.2, nrow = k),
      diag(x = 0.3, nrow = k),
      diag(x = 0.4, nrow = k),
      diag(x = 0.5, nrow = k)
    )
    psi <- list(iden)
    nu <- list(null_vec)
    lambda <- list(iden)
    theta <- list(diag(x = 0.50, nrow = k))
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
    ssm <- SimSSMIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      nu = nu,
      lambda = lambda,
      theta = theta,
      type = 0,
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
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 1
    ssm <- SimSSMIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      nu = nu,
      lambda = lambda,
      theta = theta,
      gamma_eta = gamma_eta,
      x = x,
      type = 1,
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
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 2
    ssm <- SimSSMIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      nu = nu,
      lambda = lambda,
      theta = theta,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = 2,
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
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)

    # Error
    testthat::test_that(
      paste(text, "error"),
      {
        testthat::expect_error(
          SimSSMIVary(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            alpha = alpha,
            beta = beta,
            psi = psi,
            nu = nu,
            lambda = lambda,
            theta = theta,
            gamma_y = gamma_y,
            gamma_eta = gamma_eta,
            x = x,
            type = 3,
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
          SimSSMIVary(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            alpha = alpha,
            beta = beta,
            psi = psi,
            nu = nu,
            lambda = lambda,
            theta = theta,
            type = 1,
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
          SimSSMIVary(
            n = n,
            mu0 = mu0,
            sigma0 = sigma0,
            alpha = alpha,
            beta = beta,
            psi = psi,
            nu = nu,
            lambda = lambda,
            theta = theta,
            type = 2,
            time = time,
            burn_in = burn_in
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-i-vary"
)
