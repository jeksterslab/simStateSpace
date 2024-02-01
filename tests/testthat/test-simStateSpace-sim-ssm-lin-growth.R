## ---- test-simStateSpace-sim-ssm-lin-growth
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    ## number of individuals
    n <- 5
    ## time points
    time <- 50
    ## dynamic structure
    p <- 2
    mu0 <- c(0.615, 1.006)
    sigma0 <- matrix(
      data = c(
        1.932,
        0.618,
        0.618,
        0.587
      ),
      nrow = p
    )
    sigma0_l <- t(chol(sigma0))
    ## measurement model
    k <- 1
    theta <- 0.50
    theta_l <- sqrt(theta)
    ## covariates
    j <- 2
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        return(
          matrix(
            data = rnorm(n = j * time),
            nrow = j
          )
        )
      }
    )
    gamma_eta <- diag(x = 0.10, nrow = p, ncol = j)
    gamma_y <- diag(x = 0.10, nrow = k, ncol = j)

    # Type 0
    ssm <- SimSSMLinGrowth(
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
    as.matrix.simstatespace(ssm, eta = TRUE)
    as.matrix.simstatespace(ssm, eta = FALSE)
    as.matrix.simstatespace(ssm, eta = TRUE, long = FALSE)
    as.matrix.simstatespace(ssm, eta = FALSE, long = FALSE)
    print.simstatespace(ssm)
    plot.simstatespace(ssm, id = 1:3, time = 0:4)
    plot.simstatespace(ssm, eta = TRUE)

    # Type 1
    ssm <- SimSSMLinGrowth(
      n = n,
      time = time,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      theta_l = theta_l,
      type = 1,
      x = x,
      gamma_eta = gamma_eta
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
    ssm <- SimSSMLinGrowth(
      n = n,
      time = time,
      mu0 = mu0,
      sigma0_l = sigma0_l,
      theta_l = theta_l,
      type = 2,
      x = x,
      gamma_eta = gamma_eta,
      gamma_y = gamma_y
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
  },
  text = "test-simStateSpace-sim-ssm-lin-growth"
)
