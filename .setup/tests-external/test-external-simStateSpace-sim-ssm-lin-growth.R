## ---- test-external-simStateSpace-sim-ssm-lin-growth
lapply(
  X = 1,
  FUN = function(i,
                 tol,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)
    n <- 100000
    mu0 <- c(0.615, 1.006)
    sigma0 <- matrix(
      data = c(
        1.932,
        0.618,
        0.618,
        0.587
      ),
      nrow = 2
    )
    theta <- 0.6
    time <- 5
    gamma_y <- matrix(data = 0.10, nrow = 1, ncol = 2)
    gamma_eta <- matrix(data = 0.10, nrow = 2, ncol = 2)
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

    data <- simStateSpace::Sim2Matrix(ssm, eta = FALSE, long = FALSE)

    model <- "
      # factor loadings
      eta0 =~ 1 * y_0 + 1 * y_1 + 1 * y_2 + 1 * y_3 + 1 * y_4
      eta1 =~ 0 * y_0 + 1 * y_1 + 2 * y_2 + 3 * y_3 + 4 * y_4
      # means of latent variables
      eta0 ~ mu0 * 1
      eta1 ~ mu1 * 1
      # variances and covariances of latent variables
      eta0 ~~ sigma00 * eta0
      eta0 ~~ sigma01 * eta1
      eta1 ~~ sigma11 * eta1
      # constrain error variance theta to be equal
      y_0 ~~ theta * y_0
      y_1 ~~ theta * y_1
      y_2 ~~ theta * y_2
      y_3 ~~ theta * y_3
      y_4 ~~ theta * y_4
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
  text = "test-external-simStateSpace-sim-ssm-lin-growth"
)
