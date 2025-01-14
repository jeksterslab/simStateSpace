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
  },
  text = "test-simStateSpace-tests"
)
