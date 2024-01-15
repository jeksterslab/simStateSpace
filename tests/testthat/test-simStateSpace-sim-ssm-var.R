## ---- test-simStateSpace-sim-ssm-var
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    k <- 3
    iden <- diag(k)
    null_vec <- rep(x = 0, times = k)
    mu0 <- null_vec
    sigma0 <- iden
    alpha <- null_vec
    beta <- diag(x = 0.5, nrow = k)
    psi <- iden
    time <- 50
    burn_in <- 10
    gamma_eta <- 0.10 * diag(k)
    x <- matrix(
      data = rnorm(n = k * (time + burn_in)),
      ncol = k
    )

    # No covariates
    ssm <- SimSSMVAR(
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
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

    # With covariates
    ssm <- SimSSMVAR(
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      gamma_eta = gamma_eta,
      x = x,
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

    # coverage - AR
    set.seed(42)
    k <- 1
    iden <- diag(k)
    null_vec <- rep(x = 0, times = k)
    mu0 <- null_vec
    sigma0 <- iden
    alpha <- null_vec
    beta <- diag(x = 0.5, nrow = k)
    psi <- iden
    time <- 50
    burn_in <- 10

    ssm <- SimSSMVAR(
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
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
  },
  text = "test-simStateSpace-sim-ssm-var"
)
