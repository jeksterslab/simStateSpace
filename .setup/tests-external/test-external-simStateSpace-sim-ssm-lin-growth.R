## ---- test-external-simStateSpace-sim-ssm-lin-growth
lapply(
  X = 1,
  FUN = function(i,
                 tol,
                 n,
                 text) {
    message(text)
    # prepare parameters
    p <- 2
    mu0 <- round(
      runif(n = p),
      digits = 3
    )
    sigma0 <- round(
      crossprod(
        matrix(
          data = runif(n = p^2) * 2 - 1,
          ncol = p
        )
      ),
      digits = 3
    )
    theta <- round(
      exp(
        x = rnorm(n = 1)
      ),
      digits = 3
    )
    time <- 5
    gamma_y <- round(
      matrix(
        data = rnorm(n = 2),
        nrow = 1,
        ncol = 2
      ),
      digits = 3
    )
    gamma_eta <- round(
      matrix(
        data = rnorm(n = 4),
        nrow = 2,
        ncol = 2
      ),
      digits = 3
    )
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        return(
          matrix(
            data = rnorm(n = 2 * time),
            ncol = 2
          )
        )
      }
    )

    # Type 0
    ssm <- simStateSpace::SimSSMLinGrowth(
      n = n,
      mu0 = mu0,
      sigma0 = sigma0,
      theta = theta,
      type = 0,
      time = time
    )

    data <- simStateSpace:::.Wide(ssm, eta = FALSE)

    model <- "
      # factor loadings
      eta0 =~ 1 * y1_0 + 1 * y1_1 + 1 * y1_2 + 1 * y1_3 + 1 * y1_4
      eta1 =~ 0 * y1_0 + 1 * y1_1 + 2 * y1_2 + 3 * y1_3 + 4 * y1_4
      # means of latent variables
      eta0 ~ mu0 * 1
      eta1 ~ mu1 * 1
      # variances and covariances of latent variables
      eta0 ~~ sigma00 * eta0
      eta0 ~~ sigma01 * eta1
      eta1 ~~ sigma11 * eta1
      # constrain error variance theta to be equal
      y1_0 ~~ theta * y1_0
      y1_1 ~~ theta * y1_1
      y1_2 ~~ theta * y1_2
      y1_3 ~~ theta * y1_3
      y1_4 ~~ theta * y1_4
    "

    fit <- lavaan::growth(
      model = model,
      data = as.data.frame(data)
    )

    estimates <- lavaan::coef(fit)
    mu0_hat <- round(
      c(estimates["mu0"], estimates["mu1"]),
      digits = 3
    )
    sigma0_hat <- round(
      matrix(
        data = c(
          estimates["sigma00"],
          estimates["sigma01"],
          estimates["sigma01"],
          estimates["sigma11"]
        ),
        nrow = 2
      ),
      digits = 3
    )

    theta_hat <- round(
      estimates["theta"],
      digits = 3
    )

    testthat::test_that(
      paste(text, "type 0 - mu0"),
      {
        testthat::expect_true(
          all(
            abs(
              mu0 - mu0_hat
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "type 0 - sigma0"),
      {
        testthat::expect_true(
          all(
            abs(
              sigma0 - sigma0_hat
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "type 0 - theta"),
      {
        testthat::expect_true(
          all(
            abs(
              theta - theta_hat
            ) <= tol
          )
        )
      }
    )
  },
  tol = 0.01,
  n = 100000,
  text = "test-external-simStateSpace-sim-ssm-lin-growth"
)
