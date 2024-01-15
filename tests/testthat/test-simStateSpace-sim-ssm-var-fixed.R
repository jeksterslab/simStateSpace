## ---- test-simStateSpace-sim-ssm-var-fixed
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
    n <- 5
    mu0 <- null_vec
    sigma0 <- iden
    alpha <- null_vec
    beta <- diag(x = 0.5, nrow = k)
    psi <- iden
    time <- 50
    burn_in <- 10
    gamma_eta <- 0.10 * diag(k)
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
    ssm <- SimSSMVARFixed(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      alpha = alpha,
      beta = beta,
      psi = psi,
      time = time,
      burn_in = burn_in
    )

    simStateSpace::as.data.frame.simstatespace(ssm, eta = TRUE)
    simStateSpace::as.data.frame.simstatespace(ssm, eta = FALSE)
    simStateSpace::as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    simStateSpace::as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = TRUE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = FALSE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    simStateSpace::print.simstatespace(ssm)
    simStateSpace::plot.simstatespace(ssm, id = 1:3, time = 1:10)
    simStateSpace::plot.simstatespace(ssm, eta = TRUE)

    # With covariates
    ssm <- SimSSMVARFixed(
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

    simStateSpace::as.data.frame.simstatespace(ssm, eta = TRUE)
    simStateSpace::as.data.frame.simstatespace(ssm, eta = FALSE)
    simStateSpace::as.data.frame.simstatespace(ssm, eta = TRUE, long = FALSE)
    simStateSpace::as.data.frame.simstatespace(ssm, eta = FALSE, long = FALSE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = TRUE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = FALSE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    simStateSpace::as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    simStateSpace::print.simstatespace(ssm)
    simStateSpace::plot.simstatespace(ssm, id = 1:3, time = 1:10)
    simStateSpace::plot.simstatespace(ssm, eta = TRUE)
  },
  text = "test-simStateSpace-sim-ssm-var-fixed"
)
