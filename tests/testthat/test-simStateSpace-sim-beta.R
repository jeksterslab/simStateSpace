## ---- test-simStateSpace-sim-beta
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    beta <- matrix(
      data = c(
        0.7, 0.5, -0.1,
        0.0, 0.6, 0.4,
        0, 0, 0.5
      ),
      nrow = 3
    )
    n <- 100000
    vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
    testthat::test_that(
      text,
      {
        testthat::expect_true(
          all(
            beta - (1 / n) * Reduce(
              f = `+`,
              x = SimBeta(
                n = n,
                beta = beta,
                vcov_beta_vec_l = vcov_beta_vec_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-beta"
)
