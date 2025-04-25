## ---- test-simStateSpace-sim-cov-diag-n
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
        sigma <- diag(3)
        vcov_sigma_diag_l <- t(chol(0.001 * diag(3)))
        testthat::expect_true(
          all(
            sigma - (1 / n) * Reduce(
              f = `+`,
              x = SimCovDiagN(
                n = n,
                sigma_diag = diag(sigma),
                vcov_sigma_diag_l = vcov_sigma_diag_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-cov-diag-n",
  n = 100000
)
