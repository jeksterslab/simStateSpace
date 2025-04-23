## ---- test-simStateSpace-sim-iota-n
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
        iota <- c(0, 0)
        vcov_iota_l <- t(chol(0.001 * diag(2)))
        testthat::expect_true(
          all(
            iota - (1 / n) * Reduce(
              f = `+`,
              x = SimIotaN(
                n = n,
                iota = iota,
                vcov_iota_l = vcov_iota_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-iota-n",
  n = 100000
)
