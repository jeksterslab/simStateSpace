## ---- test-simStateSpace-sim-beta-n-2
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 n) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        beta <- matrix(
          data = c(
            0.7, 0.5, -0.1,
            0.0, 0.6, 0.4,
            0, 0, 0.5
          ),
          nrow = 3
        )
        vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
        testthat::expect_true(
          all(
            beta - (1 / n) * Reduce(
              f = `+`,
              x = SimBetaN2(
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
  text = "test-simStateSpace-sim-beta-n-2",
  n = 100000
)
