## ---- test-simStateSpace-sim-mvn
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
        sigma_l <- t(chol(0.001 * diag(2)))
        testthat::expect_true(
          all(
            mu - (1 / n) * Reduce(
              f = `+`,
              x = SimMVN(
                n = n,
                mu = mu,
                sigma_l = sigma_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-mvn",
  n = 100000
)
