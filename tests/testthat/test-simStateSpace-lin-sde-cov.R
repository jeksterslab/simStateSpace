## ---- test-simStateSpace-lin-sde-cov
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
        phi <- matrix(
          data = c(
            -0.10,
            0.05,
            0.05,
            -0.10
          ),
          nrow = 2
        )
        sigma <- matrix(
          data = c(
            2.79,
            0.06,
            0.06,
            3.27
          ),
          nrow = 2
        )
        testthat::expect_true(
          all(
            (
              c(
                LinSDECov(
                  phi = phi,
                  sigma = sigma
                )
              ) - c(
                19.2,
                10.5,
                10.5,
                21.6
              )
            ) < tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-cov",
  tol = 0.000001
)
