## ---- test-simStateSpace-sim-alpha-n
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
        alpha <- c(0, 0)
        vcov_alpha_l <- t(chol(0.001 * diag(2)))
        testthat::expect_true(
          all(
            alpha - (1 / n) * Reduce(
              f = `+`,
              x = SimAlphaN(
                n = n,
                alpha = alpha,
                vcov_alpha_l = vcov_alpha_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-alpha-n",
  n = 100000
)
