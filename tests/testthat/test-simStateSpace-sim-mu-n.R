## ---- test-simStateSpace-sim-mu-n
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
        mu <- c(0, 0)
        vcov_mu_l <- t(chol(0.001 * diag(2)))
        testthat::expect_true(
          all(
            mu - (1 / n) * Reduce(
              f = `+`,
              x = SimMuN(
                n = n,
                mu = mu,
                vcov_mu_l = vcov_mu_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-mu-n",
  n = 100000
)
