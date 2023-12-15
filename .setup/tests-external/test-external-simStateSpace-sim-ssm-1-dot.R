## ---- test-external-simStateSpace-sim-ssm-1-dot
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
    mu0 <- null_vec
    sigma0_sqrt <- iden_sqrt
    alpha <- null_vec
    beta <- diag(x = 0.50, nrow = k)
    psi_sqrt <- iden_sqrt
    nu <- null_vec
    lambda <- iden
    theta_sqrt <- chol(diag(x = 0.50, nrow = k))
    time <- 50
    burn_in <- 0
    gamma_eta <- 0.10 * diag(k)
    gamma_y <- 0 * diag(k)
    x <- matrix(
      data = rnorm(n = k * (time + burn_in)),
      ncol = k
    )

    # type = 1

    set.seed(42)
    simssm1 <- simStateSpace:::.SimSSM1(
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
      time = time,
      burn_in = burn_in
    )

    set.seed(42)
    simssm_type1 <- simStateSpace::SimSSM(
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

    testthat::test_that(
      paste(text, "type = 1"),
      {
        testthat::expect_true(
          identical(
            simssm1,
            simssm_type1
          )
        )
      }
    )

    # type = 2

    set.seed(42)
    simssm2 <- simStateSpace:::.SimSSM2(
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
      time = time,
      burn_in = burn_in
    )

    set.seed(42)
    simssm_type2 <- simStateSpace::SimSSM(
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

    testthat::test_that(
      paste(text, "type = 2"),
      {
        testthat::expect_true(
          identical(
            simssm2,
            simssm_type2
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "type 1 vs. 2"),
      {
        testthat::expect_true(
          identical(
            simssm1,
            simssm2
          )
        )
      }
    )
  },
  text = "test-external-simStateSpace-sim-ssm-1-dot"
)
