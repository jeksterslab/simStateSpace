## ---- test-simStateSpace-sim-ssm-var-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    # In this example, beta varies across individuals
    set.seed(42)
    k <- 3
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
    time <- 50
    burn_in <- 10
    gamma_eta <- list(0.10 * diag(k))
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

    # No covariates
    ssm <- SimSSMVARIVary(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
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
    plot.simstatespace(ssm, id = 1:3, time = 1:10)
    plot.simstatespace(ssm, eta = TRUE)

    # With covariates
    ssm <- SimSSMVARIVary(
      n = n,
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

    as.data.frame.simstatespace(ssm, eta = TRUE)
    as.data.frame.simstatespace(ssm, eta = FALSE)
    as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = 1:10)
    plot.simstatespace(ssm, eta = TRUE)
  },
  text = "test-simStateSpace-sim-ssm-var-i-vary"
)
