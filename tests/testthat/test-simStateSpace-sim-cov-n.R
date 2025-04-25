## ---- test-simStateSpace-sim-cov-n
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
        sigma <- matrix(
          data = c(
            1.0, 0.5, 0.3,
            0.5, 1.0, 0.4,
            0.3, 0.4, 1.0
          ),
          nrow = 3
        )
        vcov_sigma_vech_l <- t(
          chol(
            0.001 * diag(3 * (3 + 1) / 2)
          )
        )
        testthat::expect_true(
          all(
            sigma - (1 / n) * Reduce(
              f = `+`,
              x = SimCovN(
                n = n,
                sigma = sigma,
                vcov_sigma_vech_l = vcov_sigma_vech_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-cov-n",
  n = 100000
)
