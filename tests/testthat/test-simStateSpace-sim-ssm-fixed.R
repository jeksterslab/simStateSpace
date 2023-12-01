## ---- test-simStateSpace-sim-ssm-fixed
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    k <- p <- 3
    iden <- diag(k)
    iden_sqrt <- chol(iden)
    null_vec <- rep(x = 0, times = k)
    n <- 5
    mu0 <- null_vec
    sigma0_sqrt <- iden_sqrt
    alpha <- null_vec
    beta <- diag(x = 0.50, nrow = k)
    psi_sqrt <- iden_sqrt
    nu <- null_vec
    lambda <- iden
    theta_sqrt <- chol(diag(x = 0.50, nrow = k))
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
    ssm <- SimSSMFixed(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      type = 0,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Type 1
    ssm <- SimSSMFixed(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      gamma_eta = gamma_eta,
      x = x,
      type = 1,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Type 2
    ssm <- SimSSMFixed(
      n = n,
      mu0 = mu0,
      sigma0_sqrt = sigma0_sqrt,
      alpha = alpha,
      beta = beta,
      psi_sqrt = psi_sqrt,
      nu = nu,
      lambda = lambda,
      theta_sqrt = theta_sqrt,
      gamma_y = gamma_y,
      gamma_eta = gamma_eta,
      x = x,
      type = 2,
      time = time,
      burn_in = burn_in
    )

    Sim2Matrix(ssm, eta = TRUE)
    Sim2Matrix(ssm, eta = FALSE)
    Sim2Matrix(ssm, eta = TRUE, long = FALSE)
    Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    # Error
    testthat::test_that(
      paste(text, "error"),
      {
        testthat::expect_error(
          SimSSMFixed(
            n = n,
            mu0 = mu0,
            sigma0_sqrt = sigma0_sqrt,
            alpha = alpha,
            beta = beta,
            psi_sqrt = psi_sqrt,
            nu = nu,
            lambda = lambda,
            theta_sqrt = theta_sqrt,
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
  },
  text = "test-simStateSpace-sim-ssm-fixed"
)
