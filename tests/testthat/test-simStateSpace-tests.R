## ---- test-simStateSpace-tests
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste(text, "stationary"),
      {
        testthat::skip_on_cran()
        testthat::expect_true(
          TestStationarity(
            matrix(
              data = c(0.5, 0.3, 0.2, 0.4),
              nrow = 2
            )
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "nonstationary"),
      {
        testthat::skip_on_cran()
        testthat::expect_false(
          TestStationarity(
            matrix(
              data = c(0.9, -0.5, 0.8, 0.7),
              nrow = 2
            )
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "stable"),
      {
        testthat::skip_on_cran()
        testthat::expect_true(
          TestStability(
            matrix(
              data = c(
                -0.357, 0.771, -0.450,
                0.0, -0.511, 0.729,
                0, 0, -0.693
              ),
              nrow = 3
            )
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "Hurwitz"),
      {
        phi <- matrix(
          data = c(
            0.10, -0.40,
            0.50, 0.20
          ),
          nrow = 2
        )
        testthat::skip_on_cran()
        testthat::expect_true(
          !TestPhiHurwitz(phi = phi)
        )
        phi <- matrix(
          data = c(
            -0.50, -0.20,
            1.00, -0.30
          ),
          nrow = 2
        )
        testthat::expect_true(
          TestPhiHurwitz(phi = phi)
        )
        testthat::expect_true(
          TestPhiHurwitz(phi = phi, eps = 1e-12)
        )
      }
    )
  },
  text = "test-simStateSpace-tests"
)
