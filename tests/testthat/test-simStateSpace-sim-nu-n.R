## ---- test-simStateSpace-sim-nu-n
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
        nu <- c(0, 0)
        vcov_nu_l <- t(chol(0.001 * diag(2)))
        testthat::expect_true(
          all(
            nu - (1 / n) * Reduce(
              f = `+`,
              x = SimNuN(
                n = n,
                nu = nu,
                vcov_nu_l = vcov_nu_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-nu-n",
  n = 100000
)
