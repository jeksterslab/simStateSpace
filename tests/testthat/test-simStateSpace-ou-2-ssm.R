## ---- test-simStateSpace-ou-2-ssm
lapply(
  X = 1,
  FUN = function(i,
                 tol,
                 text) {
    message(text)
    # prepare parameters
    set.seed(42)

    p <- k <- 2
    mu <- c(5.76, 5.18)
    phi <- matrix(
      data = c(0.10, -0.05, -0.05, 0.10),
      nrow = p
    )
    sigma <- matrix(
      data = c(2.79, 0.06, 0.06, 3.27),
      nrow = p
    )
    delta_t <- 0.10

    ssm <- OU2SSM(
      mu = mu,
      phi = phi,
      sigma = sigma,
      delta_t = delta_t
    )

    alpha <- c(0.03159928, 0.02296420)
    beta <- matrix(
      data = c(0.99006221, 0.00495027, 0.00495027, 0.99006221),
      nrow = p
    )
    psi <- matrix(
      data = c(0.27626309, 0.00743546, 0.00743546, 0.32378627),
      nrow = p
    )

    testthat::test_that(
      paste(text, "alpha"),
      {
        testthat::expect_true(
          all(
            abs(
              as.vector(
                ssm$alpha - alpha
              )
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "beta"),
      {
        testthat::expect_true(
          all(
            abs(
              as.vector(
                ssm$beta - beta
              )
            ) <= tol
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "psi"),
      {
        testthat::expect_true(
          all(
            abs(
              as.vector(
                ssm$psi - psi
              )
            ) <= tol
          )
        )
      }
    )
  },
  tol = 0.001,
  text = "test-simStateSpace-ou-2-ssm"
)
