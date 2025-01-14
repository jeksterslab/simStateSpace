## ---- test-simStateSpace-sim-ssm-var-i-vary
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # prepare parameters
        # In this example, beta varies across individuals.
        set.seed(42)
        ## number of individuals
        n <- 5
        ## time points
        time <- 50
        ## dynamic structure
        p <- 3
        mu0 <- list(
          rep(x = 0, times = p)
        )
        sigma0 <- 0.001 * diag(p)
        sigma0_l <- list(
          t(chol(sigma0))
        )
        alpha <- list(
          rep(x = 0, times = p)
        )
        beta <- SimBetaN(
          n = 5,
          beta = 0.1 * diag(p),
          vcov_beta_vec_l = t(
            chol(
              0.1 * diag(p^2)
            )
          )
        )
        psi <- 0.001 * diag(p)
        psi_l <- list(
          t(chol(psi))
        )
        ## covariates
        j <- 2
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            matrix(
              data = stats::rnorm(n = time * j),
              nrow = j,
              ncol = time
            )
          }
        )
        gamma <- list(
          diag(x = 0.10, nrow = p, ncol = j)
        )

        # Type 0
        ssm <- SimSSMVARIVary(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          alpha = alpha,
          beta = beta,
          psi_l = psi_l,
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
        ssm <- SimSSMVARIVary(
          n = n,
          time = time,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          alpha = alpha,
          beta = beta,
          psi_l = psi_l,
          type = 1,
          x = x,
          gamma = gamma
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

        testthat::expect_true(
          TRUE
        )
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-var-i-vary"
)
