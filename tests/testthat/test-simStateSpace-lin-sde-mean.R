## ---- test-simStateSpace-lin-sde-mean
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        iota <- c(0.317, 0.230)
        phi <- matrix(
          data = c(
            -0.10,
            0.05,
            0.05,
            -0.10
          ),
          nrow = 2
        )
        testthat::expect_true(
          all(
            (
              c(
                LinSDEMean(
                  phi = phi,
                  iota = iota
                )
              ) - c(
                5.76,
                5.18
              )
            ) < tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-mean",
  tol = 0.000001
)
